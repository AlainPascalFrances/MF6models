"""
Download MOD16A2GF (gap-filled 8-day ET, 500 m, Collection 6.1) over the CdL
watershed, for comparison with the model's UZF actual ET (calibration obs group).

AUTH — the Earthdata secret NEVER passes through the assistant.  Provide it via
EITHER (pick one):
  (A) a token/netrc file  C:\\Users\\<you>\\_netrc :
          machine urs.earthdata.nasa.gov
          login    <earthdata-username>
          password <earthdata-password-or-token>
  (B) environment variables set in your shell before running:
          set EARTHDATA_USERNAME=<user>
          set EARTHDATA_PASSWORD=<pass>
earthaccess.login() auto-detects both.

Run (in the flopy env):
  conda run -p C:/miniconda3/envs/flopy python download_mod16.py
"""
import earthaccess

OUT_DIR   = r"Y:\RS\MODIS4CDL"
SHORT     = "MOD16A2GF"
VERSION   = "061"
BBOX      = (-8.8665, 38.7683, -8.7771, 38.8715)   # W, S, E, N  (CdL watershed, padded)
TEMPORAL  = ("2000-01-01", "2026-12-31")

def main():
    auth = earthaccess.login()                     # reads _netrc / env vars; no interactive prompt needed
    if not getattr(auth, "authenticated", True):
        raise SystemExit("Earthdata authentication failed — check your _netrc / env vars.")
    print(f">> Searching {SHORT} v{VERSION} over {BBOX} for {TEMPORAL} …")
    results = earthaccess.search_data(
        short_name=SHORT, version=VERSION,
        bounding_box=BBOX, temporal=TEMPORAL,
    )
    print(f"   {len(results)} granule(s) found.")
    if not results:
        return
    files = earthaccess.download(results, OUT_DIR)
    print(f">> Downloaded {len(files)} file(s) to {OUT_DIR}")

if __name__ == "__main__":
    main()
