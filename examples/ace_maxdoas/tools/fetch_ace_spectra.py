#!/usr/bin/env python3
"""
Recreates examples/ace_maxdoas/data/ from the public ACE MAX-DOAS data set.

    Benavent, N., Garcia-Nieto, D., Cuevas, C. A. and Saiz-Lopez, A. (2020). Raw spectra
    measurements of scattered sunlight collected using a MAX-DOAS instrument in the austral
    summer of 2016/17 during the Antarctic Circumnavigation Expedition (ACE). Version 1.0,
    Zenodo, https://doi.org/10.5281/zenodo.3827443. License: CC BY 4.0.

The spectra are inside multi-gigabyte zip files; this script reads only the entries it
needs, with HTTP range requests (about 15 MB in all). Python 3 standard library only.

For 7 April 2017 (day 97, the ship was in the Atlantic off Portugal) and the UV
spectrometer (central wavelength 358 nm) it takes one complete elevation scan per hour,
07:00-18:00 UTC, and the zenith spectrum closest to local noon, and converts them:

  * the 19 fibre columns of each spectrum are summed into one spectrum;
  * the spectra are written in the MFC STD text format that QDOAS reads, with the
    elevation angle, the ship's position (from the data set's GPS files), the number of
    scans and the exposure time in the header;
  * the per-pixel wavelengths of the instrument's own calibration (Hg-Ne lines, identical
    in every file) are written once, to calibration/ace_uv_instrument.clb.

Nothing else is changed: the counts are the published counts (offset and dark current
were already subtracted by the instrument software).

    python3 fetch_ace_spectra.py [output_dir]      # default: ../data
"""
import datetime, io, math, os, sys, urllib.request, zipfile

RECORD = "https://zenodo.org/api/records/3827443/files/{}/content"
ZIPNAME = "ace_maxdoas_spectra-2017-04.zip"
DAY, BAND = 97, "358"
DATE = datetime.date(2017, 1, 1) + datetime.timedelta(days=DAY - 1)      # 2017-04-07
HOURS = range(7, 19)


class HttpFile(io.RawIOBase):
    """A read-only, seekable file over HTTP range requests, with a small block cache."""

    def __init__(self, url, block=1 << 20):
        self.url, self.pos, self.block, self.cache = url, 0, block, {}
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD")) as r:
            self.size = int(r.headers["Content-Length"])

    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def _block(self, k):
        if k not in self.cache:
            lo = k * self.block
            hi = min(self.size, lo + self.block) - 1
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={lo}-{hi}"})
            with urllib.request.urlopen(req) as r:
                self.cache[k] = r.read()
            if len(self.cache) > 64:
                self.cache.pop(next(iter(self.cache)))
        return self.cache[k]

    def read(self, n=-1):
        n = self.size - self.pos if n < 0 else min(n, self.size - self.pos)
        out = bytearray()
        while n > 0:
            k, o = divmod(self.pos, self.block)
            b = self._block(k)[o:o + n]
            out += b
            self.pos += len(b)
            n -= len(b)
        return bytes(out)

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def hms(v):
    v = int(v)
    return v // 10000, (v // 100) % 100, v % 100


def degrees(dm):
    """The GPS files give positions as degrees.minutes (40.01452 is 40 deg 01.452 min)."""
    d = int(abs(dm))
    return math.copysign(d + (abs(dm) - d) * 100 / 60, dm)


def gps_track():
    with urllib.request.urlopen(RECORD.format("ace_maxdoas_gps.zip")) as r:
        z = zipfile.ZipFile(io.BytesIO(r.read()))
    name = next(n for n in z.namelist() if n.endswith(f"GPS_J{DAY}.txt"))
    track = []
    for line in z.read(name).decode().splitlines()[1:]:
        f = line.strip().rstrip(",").split(",")
        if len(f) >= 3:
            h, m, s = hms(f[0])
            track.append((h * 3600 + m * 60 + s, degrees(float(f[1])), degrees(float(f[2]))))
    return sorted(track)


def position(track, t):
    return min(track, key=lambda p: abs(p[0] - t))[1:]


def main(out):
    print(f"reading the directory of {ZIPNAME} (range requests)...")
    z = zipfile.ZipFile(HttpFile(RECORD.format(ZIPNAME), block=4 << 20))
    base = f"ace_maxdoas_spectra-2017-04/MAXDOAS/{BAND}/J{DAY}/"
    names = set(z.namelist())

    # LiveInfo_*.358: one file per elevation scan, one row per elevation angle
    scans = []
    for n in sorted(x for x in names if x.startswith(base + "LiveInfo_")):
        rows = z.read(n).decode().strip().splitlines()
        hdr = rows[0].rstrip(",").split(",")
        recs = [dict(zip(hdr, r.rstrip(",").split(","))) for r in rows[1:]]
        scans.append(recs)
    print(f"{len(scans)} elevation scans on day {DAY}")

    def start(recs):
        h, m, s = hms(recs[0]["StartTime"])
        return h * 3600 + m * 60 + s

    chosen = []
    for hour in HOURS:
        cands = [s for s in scans if len(s) >= 12 and start(s) >= hour * 3600]
        if cands:
            chosen.append(min(cands, key=start))
    # the zenith spectrum closest to local noon (13 W: about 12:52 UTC) is the reference
    zen = [r for s in scans for r in s if float(r["Elevation"]) == 90]
    noon = min(zen, key=lambda r: float(r["EndSZA"]))

    track = gps_track()
    os.makedirs(os.path.join(out, "spectra"), exist_ok=True)
    os.makedirs(os.path.join(out, "reference"), exist_ok=True)
    os.makedirs(os.path.join(out, "calibration"), exist_ok=True)
    lam_written = None

    def convert(rec, folder):
        nonlocal lam_written
        elev = int(float(rec["Elevation"]))
        entry = f"{base}Atmos/{rec['Filename']}_{elev:02d}.{BAND}" if elev >= 0 else f"{base}Atmos/{rec['Filename']}_{elev}.{BAND}"
        if entry not in names:
            raise SystemExit(f"missing {entry}")
        lam, counts, nscans, expo = [], [], 0, 0.0
        for line in z.read(entry).decode().splitlines():
            f = line.strip().rstrip(",").split(",")
            if len(f) < 24:
                continue
            lam.append(float(f[1]))
            nscans, expo = int(float(f[3])), float(f[4])
            counts.append(sum(int(float(x)) for x in f[5:24]))
        if lam_written is None:
            with open(os.path.join(out, "calibration", "ace_uv_instrument.clb"), "w") as io_:
                io_.writelines(f"{x:.6f}\n" for x in lam)
            lam_written = lam
        elif max(abs(a - b) for a, b in zip(lam, lam_written)) > 1e-6:
            raise SystemExit(f"{entry}: different wavelength calibration")
        t1 = start([rec])
        h2, m2, s2 = hms(rec["EndTime"])
        t2 = h2 * 3600 + m2 * 60 + s2
        lat, lon = position(track, (t1 + t2) // 2)
        sza = 0.5 * (float(rec["StartSZA"]) + float(rec["EndSZA"]))
        tag = f"{DATE:%Y%m%d}_{h2:02d}{m2:02d}{s2:02d}_e{elev:+03d}"
        path = os.path.join(out, folder, tag + ".STD")
        with open(path, "w", newline="\n") as io_:
            io_.write("GDBGMNUP\n1\n%d\n" % len(counts))
            io_.writelines(f"{c}\n" for c in counts)
            io_.write(f"elev {elev:+d} SZA {sza:.1f}\n")                          # spectrum name
            io_.write("Princeton SP500i + PIXIS 400B, UV\n")
            io_.write("ACE MAX-DOAS, R/V Akademik Tryoshnikov\n")
            io_.write(f"{DATE:%m/%d/%Y}\n")
            io_.write("%02d:%02d:%02d\n" % (t1 // 3600, (t1 // 60) % 60, t1 % 60))
            io_.write("%02d:%02d:%02d\n" % (h2, m2, s2))
            io_.write("0\n0\n")
            io_.write(f"SCANS {nscans}\n")
            io_.write(f"INT_TIME {expo:.4f}\n")
            io_.write("SITE ACE cruise\n")
            io_.write(f"LONGITUDE {lon:.5f}\nLATITUDE {lat:.5f}\n")
            io_.write(f"ElevationAngle = {elev}\n")
            io_.write(f"ExposureTime = {1000 * expo:.1f}\n")
            io_.write(f"NumScans = {nscans}\n")
            io_.write(f"Temperature = {rec['Temperature']}\n")
            io_.write(f"SZA = {sza:.3f}\n")
            io_.write(f"Source = doi:10.5281/zenodo.3827443 MAXDOAS/{BAND}/J{DAY}/Atmos/{entry.split('/')[-1]} (CC BY 4.0)\n")
        return path

    n = 0
    for s in chosen:
        for rec in s:
            convert(rec, "spectra")
            n += 1
    ref = convert(noon, "reference")
    print(f"wrote {n} spectra, reference {os.path.basename(ref)}, into {out}")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "data")))
