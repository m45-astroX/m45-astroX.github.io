'''
    Program : calculate_star_temperature.py
    Author  : Yuma Aoki (Kindai Univ.)
    Date    : 2026.05.04
    Version : 1.0
'''

from pathlib import Path
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u


### 1. 設定 ###
GRID = "ck04models"

# フィット範囲
FIT_MIN = 2000.0
FIT_MAX = 3000.0

# 温度関連設定
T_LOW = 3000.0
T_HIGH = 20000.0
N_TEMP = 100

# 恒星環境設定
METALLICITY = 0.0
LOGG = 4.0
EBV = 0.0

# スペクトル関連設定
BIN_WIDTH = 25.0
SNR_MIN = 3.0
Y_PERCENTILE = 99

# モデル関連設定
TRDS_DIR = Path("/content/trds")
OBS_PATTERN = "*aspec*.fits"


### 2. リファレンスデータの設定 ###
catalog_files = sorted(TRDS_DIR.rglob(f"grid/{GRID}/catalog.fits"))

if not catalog_files:
    raise FileNotFoundError(
        f"{GRID} の catalog.fits が見つかりません。"
    )

catalog_file = catalog_files[0]
cdbs_root = catalog_file.parent.parent.parent

os.environ["PYSYN_CDBS"] = str(cdbs_root)

### 3. HST のデータ解析に必要なパッケージをインポート ###
import stsynphot as stsyn
from stsynphot.config import overwrite_synphot_config
from stsynphot.catalog import grid_to_spec
overwrite_synphot_config(str(cdbs_root))


### 4. 観測データのインポート ###
obs_files = sorted(Path(".").glob(OBS_PATTERN))

if not obs_files:
    raise FileNotFoundError(f"{OBS_PATTERN} が見つかりません。")

obs_file = obs_files[0]
#print(f"Reading observation: {obs_file}")

with fits.open(obs_file) as hdul:
    data = hdul[1].data

    wave = np.ravel(data["WAVELENGTH"])
    flux = np.ravel(data["FLUX"])

    if "SNR" in data.columns.names:
        snr = np.ravel(data["SNR"])
    else:
        snr = np.full_like(flux, np.inf)


### 5. データのスクリーニング ###
mask = np.isfinite(wave)
mask &= np.isfinite(flux)
mask &= np.isfinite(snr)
mask &= flux > 0
mask &= snr > SNR_MIN

wave = wave[mask]
flux = flux[mask]

order = np.argsort(wave)
wave = wave[order]
flux = flux[order]

#print(f"Observed wavelength range: {wave.min():.1f} - {wave.max():.1f} Å")


### 6. スペクトルのビニング ###
wave_min = np.floor(wave.min() / BIN_WIDTH) * BIN_WIDTH
wave_max = np.ceil(wave.max() / BIN_WIDTH) * BIN_WIDTH

edges = np.arange(wave_min, wave_max + BIN_WIDTH, BIN_WIDTH)

bin_wave = []
bin_flux = []

for i in range(len(edges) - 1):
    lo = edges[i]
    hi = edges[i + 1]

    in_bin = (wave >= lo) & (wave < hi)

    if np.any(in_bin):
        bin_wave.append(0.5 * (lo + hi))
        bin_flux.append(np.nanmedian(flux[in_bin]))

bin_wave = np.array(bin_wave)
bin_flux = np.array(bin_flux)

mask = np.isfinite(bin_flux)
mask &= bin_flux > 0

bin_wave = bin_wave[mask]
bin_flux = bin_flux[mask]


### 7. フィット範囲の指定 ###
fit_mask = bin_wave >= FIT_MIN
fit_mask &= bin_wave <= FIT_MAX
fit_mask &= bin_flux > 0

x_fit = bin_wave[fit_mask]
y_fit = bin_flux[fit_mask]

if len(x_fit) < 5:
    raise RuntimeError(
        "フィット点が少なすぎます。\n"
        f"観測データの範囲: {wave.min():.0f} - {wave.max():.0f} Å\n"
        f"指定した範囲: {FIT_MIN:.0f} - {FIT_MAX:.0f} Å"
    )


### 8. モデルの設定 ###
def get_model_flux(wave_A, temperature):
    spectrum = grid_to_spec(
        GRID,
        temperature,
        METALLICITY,
        LOGG,
    )

    if EBV != 0.0:
        spectrum = spectrum * stsyn.ebmvx("mwavg", EBV)

    wavelength = wave_A * u.AA
    model_flux = spectrum(wavelength)

    model_flux = model_flux.to(
        u.erg / u.s / u.cm**2 / u.AA,
        equivalencies=u.spectral_density(wavelength),
    )

    return np.asarray(model_flux.value)

def get_log_scale(obs_flux, model_flux):
    mask = np.isfinite(obs_flux)
    mask &= np.isfinite(model_flux)
    mask &= obs_flux > 0
    mask &= model_flux > 0

    obs_log = np.log10(obs_flux[mask])
    model_log = np.log10(model_flux[mask])

    return np.nanmedian(obs_log - model_log)

def get_score(temperature):
    try:
        model_flux = get_model_flux(x_fit, temperature)
    except Exception:
        return np.inf

    mask = np.isfinite(y_fit)
    mask &= np.isfinite(model_flux)
    mask &= y_fit > 0
    mask &= model_flux > 0

    if np.sum(mask) < 5:
        return np.inf

    log_scale = get_log_scale(y_fit[mask], model_flux[mask])

    obs_log = np.log10(y_fit[mask])
    model_log = np.log10(model_flux[mask]) + log_scale

    residual = obs_log - model_log

    return np.nanmean(residual**2)


### 9. 温度の計算 ###
temperatures = np.linspace(T_LOW, T_HIGH, N_TEMP)

scores = []

for temperature in temperatures:
    score = get_score(temperature)
    scores.append(score)

scores = np.array(scores)

best_index = np.nanargmin(scores)
best_temperature = temperatures[best_index]

model_flux = get_model_flux(bin_wave, best_temperature)

model_flux_fit_range = get_model_flux(x_fit, best_temperature)
log_scale = get_log_scale(y_fit, model_flux_fit_range)

model_flux = model_flux * 10**log_scale


### 10. 結果の出力 ###
print(f"温度 = {best_temperature:.0f} K")


### 11. エラー出力 ###
if best_temperature == T_LOW:
    print("Warning: Teff が下限に達しています。")

if best_temperature == T_HIGH:
    print("Warning: Teff が上限に達しています。")


### 12. スペクトルのプロット ###
plt.figure(figsize=(10, 5))

plt.plot(
    bin_wave,
    bin_flux,
    lw=1.2,
    label=f"Observed spectrum ({BIN_WIDTH:.0f} Å bins)",
)

plt.plot(
    bin_wave,
    model_flux,
    lw=2.0,
    label=f"{GRID}: Teff = {best_temperature:.0f} K",
)

plot_mask = bin_wave >= FIT_MIN
plot_mask &= bin_wave <= FIT_MAX
plot_mask &= bin_flux > 0

ymax = np.nanpercentile(bin_flux[plot_mask], Y_PERCENTILE)

plt.xlim(FIT_MIN, FIT_MAX)
plt.ylim(0, ymax)

plt.xlabel(r"Wavelength [$\AA$]")
plt.ylabel(r"Flux [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
