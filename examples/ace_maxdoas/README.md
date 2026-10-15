# Real-data example: MAX-DOAS on the Antarctic Circumnavigation Expedition

This example analyses **real, publicly available** MAX-DOAS spectra with **public
laboratory cross sections**, in a QDOAS project that QDOAS itself also runs.

```bash
julia --project examples/ace_maxdoas/prepare.jl        # once; needs internet (about 7 MB)
julia --project -t auto examples/run_gui.jl ace         # the GUI on this project
julia --project -t auto bin/minidoas.jl -c examples/ace_maxdoas/work/ace_maxdoas.xml -o ace_results
```

## The spectra (`data/`)

The spectra were measured by the MAX-DOAS instrument of the Instituto de Química Física
Rocasolano (CSIC) on the research vessel *Akademik Tryoshnikov* during the Antarctic
Circumnavigation Expedition (ACE), and published by

> Benavent, N., Garcia-Nieto, D., Cuevas, C. A. and Saiz-Lopez, A. (2020). *Raw spectra
> measurements of scattered sunlight collected using a MAX-DOAS (Multi-Axis Differential
> Optical Absorption Spectroscopy) instrument in the austral summer of 2016/17 during the
> Antarctic Circumnavigation Expedition (ACE).* Version 1.0 [Data set]. Zenodo.
> <https://doi.org/10.5281/zenodo.3827443>

under the **Creative Commons Attribution 4.0 International license (CC BY 4.0,
<https://creativecommons.org/licenses/by/4.0/>)**. The data set's authors are not
associated with MiniDOAS.jl and do not endorse it.

This folder holds a small part of that data set: the UV spectrometer (Princeton Instruments
SP500i with a PIXIS 400B CCD, 1340 pixels, 314.6–400.9 nm) on **7 April 2017**, when the
ship was in the Atlantic off Portugal (about 40°N, 12°W):

* `spectra/`: 144 spectra, one complete elevation scan per hour from 07:00 to 18:00 UTC,
  each scan at 12 viewing elevations (−2°, 0°, 1°, 2°, 3°, 4°, 6°, 10°, 15°, 30°, 70° and 90°);
* `reference/`: the zenith spectrum closest to local noon (12:51 UTC, solar zenith angle
  34.5°), used as the Fraunhofer reference;
* `calibration/ace_uv_instrument.clb`: the instrument's own wavelength calibration
  (from Hg–Ne lines, as published with every spectrum).

**Changes made to the published data** (as CC BY 4.0 asks us to state), all by
[`tools/fetch_ace_spectra.py`](tools/fetch_ace_spectra.py), which recreates this folder
from Zenodo:

* the 19 optical-fibre columns of each spectrum were added up into one spectrum; the counts
  are otherwise the published ones (offset and dark current had been subtracted by the
  instrument software);
* the spectra were written in the MFC STD text format that QDOAS reads, with the elevation
  angle, the solar zenith angle, the number of scans and the exposure time from the data
  set's `LiveInfo` files, and the ship's position from its GPS files (which give it in
  degrees and minutes, converted to decimal degrees); the last header line of each file
  names the original file;
* the wavelength column, identical in every file, was written once to the calibration file.

## What `prepare.jl` does (`work/`, not in the repository)

1. **Downloads** the cross sections from the MPI-Mainz UV/VIS Spectral Atlas (Keller-Rudek
   et al., 2013) and the solar spectrum of Chance and Kurucz (2010). They are not
   redistributed here; see the list below for whom to cite.
2. **Calibrates the wavelengths and the slit function** with `calibrate_wavelength`, in the
   spirit of QDOAS's "Kurucz" calibration: in eight 10 nm sub-windows from 318 to 398 nm it fits the high-resolution solar
   spectrum, convolved with a Gaussian slit, to the noon zenith spectrum (plus a polynomial,
   ozone and a Ring term), and keeps a smooth shift (2nd order in λ) and slit width (1st
   order). The instrument's calibration is off by +0.06 to +0.33 nm from the solar
   spectrum's vacuum wavelengths; the slit's full width at half maximum is 0.37 nm.
3. **Convolves each cross section** with that slit on the calibrated grid (`convolve`; air
   wavelengths converted to vacuum with `air_to_vacuum` where a data set uses air), and
   computes a **Ring spectrum** from the solar spectrum with `ring_spectrum`, QDOAS's own
   Ring algorithm (rotational Raman scattering of N₂ and O₂ at 250 K after Chance and Spurr,
   1997).
4. **Writes the QDOAS project** `work/ace_maxdoas.xml` (`write_qdoas_project`) with two
   analysis windows:

| window | interval | absorbers | polynomial, offset | shift and stretch |
|---|---|---|---|---|
| NO2 | 338–370 nm | NO₂, O₃, O₄, BrO, HCHO, Ring | order 5, linear offset of order 1 | on *Spectrum* |
| O4 | 352–384 nm | O₄, NO₂, O₃, BrO, Ring | order 5, linear offset of order 1 | on *Spectrum* |

The wavelength convention of the strong absorbers was checked on these spectra by fitting
each with a free shift of its own: O₄, O₃ and the Ring spectrum come out within 0.03 nm of
zero (the O₄ data turned out to be in vacuum wavelengths; assuming air would have shifted
them by 0.10 nm). NO₂, BrO and HCHO are too weak in these spectra to be checked this way;
their conventions are taken from the literature (vacuum, vacuum and air).

The results are **differential** slant columns relative to the noon zenith spectrum. This is
a demonstration set-up: a careful analysis would add an I₀ correction for ozone, several
ozone temperatures, a measured slit function and quality filtering (clouds, the −2° view of
the sea surface).

### What the analysis gives

Run on 15 October 2026 with the cross sections and solar spectrum as downloaded that day
(`prepare.jl` then `analyze_spectra`; 0.7 s for the 144 spectra on four threads):

* **Calibration**: the instrument's wavelengths are off by +0.05 nm at 323 nm to +0.34 nm at
  393 nm from the solar spectrum's vacuum wavelengths; the slit's FWHM is 0.36-0.37 nm.
* **Fit quality**: all 144 spectra fitted, in 2-3 iterations; RMS 4.5 × 10⁻⁴ to 2.8 × 10⁻³
  (median 6.4 × 10⁻⁴), largest for the −2° and 0° views; fitted spectrum shifts within
  ±0.025 nm (median −0.004 nm), so the calibration holds.
* **O₄** (O4 window, relative error about 2 %) grows towards the horizon, from about 0 in the
  zenith at noon (the reference) to 3.3 × 10⁴³ molecules² cm⁻⁵ at 4-6° (medians over the day),
  levelling off and dropping slightly at 0-2° and at −2° (the sea surface), as visibility limits
  the light path near the horizon. The geometric approximation at 30° elevation,
  `geometric_vcd(dSCD(30°) - dSCD(zenith), 30)`, gives an O₄ vertical column of
  **1.28-1.43 × 10⁴³ molecules² cm⁻⁵ between 10 and 15 UTC**, as expected for O₄ at sea level
  (about 1.2-1.4 × 10⁴³); the approximation fails, as it should, at large solar zenith angles.
  At 16 UTC all elevations collapse to about 1.5 × 10⁴³, as clouds or fog do.
* **NO₂** in the zenith follows the sun: −1.4 to +0.2 × 10¹⁵ molecules cm⁻² from 10 to 16 UTC,
  1.1 × 10¹⁶ at 07:03 UTC (SZA 84°) and 8.8 × 10¹⁵ at 18:02 UTC (SZA 77°), most of it in the
  stratosphere (larger in the evening at a smaller SZA, the expected diurnal increase of
  stratospheric NO₂). The off-axis views add almost nothing: tropospheric NO₂ from the 15° view
  is (−8 to +6) × 10¹⁴ molecules cm⁻², consistent with zero (its error is about 4 × 10¹⁴) - clean
  marine air.
* **O₃** reaches 3.7 × 10¹⁹ molecules cm⁻² at the largest SZA (stratospheric, NO2 window);
  **BrO** and **HCHO** are not detected above their errors; the **Ring** coefficient goes from
  about 0 (the noon zenith) to −0.05 (more Raman filling-in than in the reference), with a 3 %
  error.

<p align="center"><img src="../../docs/img/plot_fit_ace.png" alt="The NO2 window fitted to the 15° spectrum of 18:02 UTC" width="85%"></p>

<p align="center"><img src="../../docs/img/ace_o4.png" alt="O4 differential slant columns through the day, one line per viewing elevation" width="90%"></p>

*Above: the NO₂ window fitted to the 15° spectrum of 18:02 UTC (`plot_fit`). Below: O₄
differential slant columns through the day, one line per viewing elevation (both made by
`docs/make_figures.jl`).*

### Cross sections and solar spectrum (downloaded by `prepare.jl`)

| symbol | data set | cite |
|---|---|---|
| NO₂ | Vandaele et al. (1998), 294 K | J. Quant. Spectrosc. Radiat. Transfer 59, 171–184, doi:10.1016/S0022-4073(97)00168-4 |
| O₃ | Serdyuchenko et al. (2014), 223 K | Atmos. Meas. Tech. 7, 625–636, doi:10.5194/amt-7-625-2014 |
| O₄ | Thalman and Volkamer (2013), 293 K | Phys. Chem. Chem. Phys. 15, 15371–15381, doi:10.1039/C3CP50968K |
| BrO | Fleischmann et al. (2004), 223 K | J. Photochem. Photobiol. A 168, 117–132, doi:10.1016/j.jphotochem.2004.03.026 |
| HCHO | Meller and Moortgat (2000), 298 K | J. Geophys. Res. 105, 7089–7101, doi:10.1029/1999JD901074 |
| solar | Chance and Kurucz (2010) | J. Quant. Spectrosc. Radiat. Transfer 111, 1289–1295, doi:10.1016/j.jqsrt.2010.01.036 |

The cross sections come from the MPI-Mainz UV/VIS Spectral Atlas: Keller-Rudek, H.,
Moortgat, G. K., Sander, R. and Sörensen, R. (2013), Earth Syst. Sci. Data 5, 365–373,
doi:10.5194/essd-5-365-2013. The solar spectrum comes from the Harvard-Smithsonian Center
for Astrophysics (<https://www.cfa.harvard.edu/atmosphere/>).

## QDOAS and MiniDOAS.jl on these spectra

QDOAS 3.7.5 (`doas_cl`) and QDOASJulia, whose analysis MiniDOAS.jl runs unchanged, were run on all
144 spectra of this project as QDOASJulia prepared it: 288 window fits and 3,888 fitted
parameters agree to a median of **3.5 × 10⁻¹¹** of their errors (at most 3.3 × 10⁻⁷), with
identical iteration counts and RMS values equal to 2.6 × 10⁻¹² relative. QDOAS took 7.4 s for the
144 spectra (1.15 s for one); the Julia analysis 0.41 s on one core and 0.18 s on four.

On the same downloads, MiniDOAS.jl's `prepare.jl` (built on `calibrate_wavelength`, `convolve`,
`ring_spectrum` and `write_qdoas_project`) reproduces the calibration of QDOASJulia's script digit
for digit. Its files differ only in how many digits are written (wavelengths by at most
5 × 10⁻⁷ nm, cross sections by at most 3 × 10⁻⁹ relative), so the 3,888 fitted parameters of the
two preparations differ by at most 6.6 × 10⁻⁵ of their errors, with identical iteration counts.
On QDOASJulia's own project, MiniDOAS.jl and QDOASJulia agree bit for bit in all 288 window fits.

**Positions and solar zenith angles.** The ship's GPS files give the position as degrees and
minutes (40.01452 is 40° 01.452′); the spectra here hold decimal degrees. The solar zenith angle
MiniDOAS.jl computes from each spectrum's time and position (`solar_position(s)`) is on average 0.18°
(at most 0.46°) smaller than the `SZA` the instrument software recorded (the header line
`SZA = ...`, the mean of `StartSZA` and `EndSZA`); the difference does not follow a time offset,
and MiniDOAS.jl's calculation reproduces the solstice and equinox declinations and the equation of
time to better than 0.01° and 0.2 min.
