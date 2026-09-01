========================================================================
 CdL MODFLOW 6 model  -  HOW TO RUN  (operational quick-start)
========================================================================
Self-contained set of MODFLOW 6 input files for the calibrated CdL
catchment model (Voronoi DISV, 6609 cells/layer x 3 layers, 557 monthly
stress periods, 1981-2026, with ATS). Everything MF6 needs is in THIS
folder. No paths to edit - all file references are relative.

------------------------------------------------------------------------
 REQUIREMENTS
------------------------------------------------------------------------
 >>> The model runs with the MF6 EXECUTABLE ALONE. <<<
 To run the simulation you need nothing but mf6.exe and the folder with the input data.
 
 - MODFLOW 6 executable (mf6.exe), v6.7.0 or later, on your PATH
   (or copy mf6.exe into this folder).  <-- the ONLY software required
 - ~300 MB free disk (the UZF file expands to ~285 MB; outputs add more).

------------------------------------------------------------------------
CdL MODFLOW 6 model — input files
------------------------------------------------------------------------

   mfsim.nam / cdl_gwf.nam   simulation + GWF model name files
   cdl_gwf.tdis(.ats)        time discretisation (557 periods + ATS)
   cdl_gwf.ims               solver (Newton + under-relaxation)
   cdl_gwf.disv              grid
   cdl_gwf.ic                initial heads
   cdl_gwf.npf / .sto        flow + storage properties
   cdl_gwf.drn/_0.drn/_1.drn seepage / western / secondary drains
   cdl_gwf.ghb               eastern general-head boundary (underflow in)
   cdl_gwf.sfr               streams
   cdl_gwf.lak(+ lak*.tab)   ponds
   cdl_gwf.mvr               water mover (SFR<->LAK<->UZF routing)
   cdl_gwf.uzf(.gz)          unsaturated zone / recharge + ET (see note)
   cdl_gwf.oc                output control
   *.obs                     observation configurations

Parameters are the CALIBRATED (PEST++ ies) values baked into the inputs.

------------------------------------------------------------------------
 STEP 1  -  DECOMPRESS THE UZF FILE  (one time, required)
------------------------------------------------------------------------
 The UZF package is shipped gzipped (cdl_gwf.uzf.gz) because uncompressed
 it is ~285 MB. MF6 needs the plain file "cdl_gwf.uzf". Unzip it first:

   PowerShell:
     $in='cdl_gwf.uzf.gz'; $out='cdl_gwf.uzf'
     $fi=[IO.File]::OpenRead($in); $fo=[IO.File]::Create($out)
     $gz=New-Object IO.Compression.GzipStream($fi,[IO.Compression.CompressionMode]::Decompress)
     $gz.CopyTo($fo); $gz.Dispose(); $fo.Dispose(); $fi.Dispose()

   or Git Bash / 7-Zip:
     gunzip -k cdl_gwf.uzf.gz          (Git Bash; -k keeps the .gz)
     7z e cdl_gwf.uzf.gz               (7-Zip)

 After this the folder must contain cdl_gwf.uzf (NOT only the .gz).

------------------------------------------------------------------------
 STEP 2  -  RUN
------------------------------------------------------------------------
 Open a terminal IN THIS FOLDER and launch mf6 with no arguments; it
 reads mfsim.nam automatically:

     cd /d E:\MF6models\CdL\input_model      (Windows cmd)
     mf6

 A full transient run takes on the order of ~30 minutes (ATS adapts the
 time step). Watch the console: it should end with "Normal termination of simulation".

------------------------------------------------------------------------
 STEP 3  -  OUTPUTS  (written into this folder by the run)
------------------------------------------------------------------------
   mfsim.lst          simulation listing (check here first if it fails)
   cdl_gwf.lst        model listing + volumetric water budget
   cdl_gwf.hds        heads (binary)
   cdl_gwf.cbc        cell-by-cell flows (binary)
   cdl_gwf.*.bud      SFR / LAK / UZF package budgets
   cdl_gwf.grb        binary grid (needed to post-process the DISV grid)
   *.csv              observation output (heads, GHB, SFR, DRN, lake stage)

------------------------------------------------------------------------
 NOTES / TROUBLESHOOTING
------------------------------------------------------------------------
 - "cdl_gwf.uzf ... file not found": you skipped Step 1 (decompress).
 - "mf6 is not recognized": mf6.exe is not on PATH; copy it here or add
   its folder to PATH.
 - Do not rename files - the name files reference them by these exact
   names.
 - Note that the MODFLOW6 input files were generated using Python scripts
   that are not provided here.
========================================================================
