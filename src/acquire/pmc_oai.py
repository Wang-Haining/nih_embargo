"""Retrieve PMC OAI-PMH front matter (pmc_fm) for each PMCID and parse the dates and flags used here.

Fields: the pmc-live event (the day the record went live), the pmc-release event (the release date
PMC reports), the NIH Manuscript Submission event, the accepted date, and the manuscript, live,
embargo and licence flags from custom-meta. All events and custom-meta pairs are kept as JSON, and
each raw record is saved gzipped to raw/. PMC asks that jobs over 100 requests run outside 05:00 to
21:00 US Eastern on weekdays; the script waits for that window. Resumable: rows are checkpointed to
the output parquet. data/pmc_records.csv was built from the outputs of this script, run on 16 and
17 October 2026.

usage: python src/acquire/pmc_oai.py <input.csv with pmid,pmcid> <output.parquet> [--now]
"""

import gzip
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd
from lxml import etree

URL = "https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/?verb=GetRecord&metadataPrefix=pmc_fm&identifier=oai:pubmedcentral.nih.gov:{}"
INTERVAL = 0.8
RAW = Path("raw")
META = [
    "pmc-prop-manuscript",
    "pmc-prop-open-access",
    "pmc-prop-olf",
    "pmc-status-embargo",
    "pmc-status-live",
    "pmc-status-released",
    "pmc-license-ref",
    "pmc-collection-title",
]


def off_peak() -> bool:
    now = datetime.now(ZoneInfo("America/New_York"))
    return now.weekday() >= 5 or now.hour >= 21 or now.hour < 5


def local(e) -> str:
    return etree.QName(e).localname


def event_date(e):
    for d in e.iter():
        if isinstance(d.tag, str) and local(d) == "date":
            iso = d.get("iso-8601-date")
            if iso:
                return pd.Timestamp(iso[:10])
            parts = {local(x): (x.text or "").strip() for x in d if isinstance(x.tag, str)}
            if parts.get("year"):
                return pd.Timestamp(
                    int(parts["year"]), int(parts.get("month") or 1), int(parts.get("day") or 1)
                )
    return None


def pub_date(e):
    parts = {local(x): (x.text or "").strip() for x in e if isinstance(x.tag, str)}
    if not parts.get("year"):
        return None
    return pd.Timestamp(int(parts["year"]), int(parts.get("month") or 1), int(parts.get("day") or 1))


def parse(pmid, pmcid, body: bytes) -> dict:
    root = etree.fromstring(body)
    row = {"pmid": pmid, "pmcid": pmcid, "error": None}
    for e in root.iter():
        if isinstance(e.tag, str) and local(e) == "error":
            row["error"] = (e.text or "").strip()
            return row
    events, accepted, meta, epub = {}, [], {}, []
    all_events, all_meta = [], []
    for e in root.iter():
        if not isinstance(e.tag, str):
            continue
        n = local(e)
        if n == "datestamp":
            row["oai_datestamp"] = (e.text or "").strip()
        elif n == "event":
            d = event_date(e)
            all_events.append(
                [e.get("event-type"), e.get("specific-use"), str(d.date()) if d is not None else None]
            )
            events.setdefault(e.get("event-type"), []).append(d)
        elif n == "date" and e.get("date-type") == "accepted":
            accepted.append(pub_date(e))
        elif n == "pub-date" and (e.get("pub-type") == "epub" or e.get("date-type") == "pub"):
            epub.append(pub_date(e))
        elif n == "custom-meta":
            k = [x for x in e if isinstance(x.tag, str) and local(x) == "meta-name"]
            v = [x for x in e if isinstance(x.tag, str) and local(x) == "meta-value"]
            if k:
                all_meta.append([k[0].text, v[0].text if v else None])
            if k and k[0].text in META:
                meta[k[0].text] = v[0].text if v else None
    first = lambda xs: min((x for x in xs if x is not None), default=None)  # noqa: E731
    row["pmc_release"] = first(events.get("pmc-release", []))
    row["pmc_live"] = first(events.get("pmc-live", []))
    row["nihms_event"] = "nihms-submitted" in events
    row["nihms_submitted"] = first(events.get("nihms-submitted", []))
    row["accepted"] = first(accepted)
    row["epub"] = first(epub)
    row["pmc_released_under_embargo"] = first(events.get("pmc-released-under-embargo", []))
    row["events_json"] = json.dumps(all_events)
    row["meta_json"] = json.dumps(all_meta)
    for k in META:
        row[k.replace("-", "_")] = meta.get(k)
    return row


def main() -> None:
    src, out = sys.argv[1], sys.argv[2]
    todo = pd.read_csv(src, dtype=str)[["pmid", "pmcid"]].dropna().drop_duplicates("pmcid")
    try:
        done = pd.read_parquet(out)
    except FileNotFoundError:
        done = pd.DataFrame()
    rows = done.to_dict("records")
    seen = set(done.pmcid) if len(done) else set()
    todo = todo[~todo.pmcid.isin(seen)]
    print(f"{len(seen)} done, {len(todo)} to fetch", flush=True)
    last = 0.0
    for i, r in enumerate(todo.itertuples(index=False), 1):
        while len(todo) > 100 and "--now" not in sys.argv and not off_peak():
            time.sleep(300)
        for attempt in range(6):
            time.sleep(max(0.0, INTERVAL - (time.monotonic() - last)))
            try:
                req = urllib.request.Request(
                    URL.format(r.pmcid.removeprefix("PMC")),
                    headers={
                        "User-Agent": "nih-embargo/1.0 (hw56@iu.edu)",
                        "Accept-Encoding": "gzip, deflate",
                    },
                )
                with urllib.request.urlopen(req, timeout=90) as resp:
                    body = resp.read()
                    if resp.headers.get("Content-Encoding", "").lower() == "gzip":
                        body = gzip.decompress(body)
                last = time.monotonic()
                RAW.mkdir(exist_ok=True)
                (RAW / f"{r.pmcid}.xml.gz").write_bytes(gzip.compress(body))
                rows.append(parse(r.pmid, r.pmcid, body))
                break
            except urllib.error.HTTPError as exc:
                last = time.monotonic()
                if 400 <= exc.code < 500 and exc.code != 429:
                    rows.append({"pmid": r.pmid, "pmcid": r.pmcid, "error": f"HTTP {exc.code}"})
                    break
                if attempt == 5:
                    rows.append({"pmid": r.pmid, "pmcid": r.pmcid, "error": f"HTTP {exc.code}"})
                time.sleep(2**attempt)
            except Exception as exc:  # retry transient network and parse errors
                last = time.monotonic()
                if attempt == 5:
                    rows.append({"pmid": r.pmid, "pmcid": r.pmcid, "error": repr(exc)[:200]})
                time.sleep(2**attempt)
        if i % 250 == 0 or i == len(todo):
            pd.DataFrame(rows).to_parquet(out, index=False)
            print(f"{i}/{len(todo)} {datetime.now():%H:%M}", flush=True)
    pd.DataFrame(rows).to_parquet(out, index=False)


if __name__ == "__main__":
    main()
