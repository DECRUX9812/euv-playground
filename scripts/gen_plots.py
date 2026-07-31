#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import matplotlib.pyplot as plt
import numpy as np
from euv_playground.multilayer import reflectivity_vs_angle, reflectivity_vs_wavelength
from euv_playground.zernike import phase_map, zernike
from euv_playground.imaging import (abbe_image_2d, contrast_vs_pitch,
                                    dense_lines, image_grating)
from euv_playground.resist import plot_contrast_curve, shot_noise

OUT = ROOT/"plots"; OUT.mkdir(exist_ok=True)
lam=np.linspace(12.5,14.5,601); fig,ax=plt.subplots(); ax.plot(lam,reflectivity_vs_wavelength(lam)); ax.axvline(13.5,color='r',ls='--'); ax.set(xlabel='Wavelength (nm)',ylabel='Reflectivity',title='Mo/Si multilayer'); fig.tight_layout(); fig.savefig(OUT/'reflectivity_wavelength.png',dpi=160); plt.close(fig)
a=np.linspace(0,10,401); fig,ax=plt.subplots(); ax.plot(a,reflectivity_vs_angle(a)); ax.set(xlabel='Angle from normal (degrees)',ylabel='Reflectivity'); fig.tight_layout(); fig.savefig(OUT/'reflectivity_angle.png',dpi=160); plt.close(fig)
g=np.linspace(-1,1,301); xx,yy=np.meshgrid(g,g); rr=np.hypot(xx,yy); tt=np.arctan2(yy,xx); terms=[(2,0,'Defocus'),(2,2,'Astigmatism'),(2,-2,'Astigmatism 45°'),(3,1,'Coma X'),(3,-1,'Coma Y'),(3,3,'Trefoil'),(4,0,'Spherical'),(4,2,'Secondary astig.'),(4,4,'Quadrafoil')]; fig,axs=plt.subplots(3,3,figsize=(8,8));
for ax,(n,m,title) in zip(axs.flat,terms): ax.imshow(zernike(n,m,rr,tt),extent=(-1,1,-1,1),cmap='RdBu'); ax.set_title(title); ax.axis('off')
fig.tight_layout(); fig.savefig(OUT/'zernike_maps.png',dpi=160); plt.close(fig)
x=np.linspace(-1000,1000,4097); mask,clear=image_grating(x,64); _,blur=image_grating(x,64,zernike_terms=[(2,0,1.5)]); fig,ax=plt.subplots(); ax.plot(x,clear,label='No aberration'); ax.plot(x,blur,label='Defocus'); ax.set_xlim(-250,250); ax.set(xlabel='Position (nm)',ylabel='Intensity'); ax.legend(); fig.tight_layout(); fig.savefig(OUT/'imaging_demo.png',dpi=160); plt.close(fig)
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,4)); plot_contrast_curve(ax1); doses=np.linspace(5,100,200); _,p,n=shot_noise(doses,10); ax2.plot(doses,p,label='photons/pixel'); ax2.set(xlabel='Dose (mJ/cm²)',ylabel='Photons per 10 nm pixel'); axn=ax2.twinx(); axn.plot(doses,100*n,color='tab:red',label='noise'); axn.set_ylabel('Noise (%)',color='tab:red'); fig.tight_layout(); fig.savefig(OUT/'resist_and_shot_noise.png',dpi=160); plt.close(fig)

n2, dx2 = 256, 1.0; mask2 = dense_lines(n2, dx2, 64)
image_cases = [("Ideal", ()), ("Defocus", ((2, 0, 0.8),)),
               ("Coma", ((3, 1, 0.8),)), ("Astigmatism", ((2, 2, 0.8),))]
fig, axs = plt.subplots(2, 2, figsize=(8, 7))
for ax, (title, aberration) in zip(axs.flat, image_cases):
    image = abbe_image_2d(mask2, dx2, zernike_terms=aberration)
    ax.imshow(image, cmap="gray", extent=(-128, 128, -128, 128), origin="lower")
    ax.set(title=title, xlabel="x (nm)", ylabel="y (nm)")
fig.tight_layout(); fig.savefig(OUT/'imaging_2d_grid.png',dpi=160); plt.close(fig)

pitches = np.arange(14.0, 81.0, 2.0)
fig, ax = plt.subplots()
for sigma, label, cutoff in [(0.0, "σ = 0", 13.5/.33),
                             (1.0, "σ = 1", 13.5/(2*.33))]:
    curve = np.asarray(contrast_vs_pitch(256, 1.0, pitches, sigma))
    ax.plot(curve[:, 0], curve[:, 1], marker=".", label=label)
    ax.axvline(cutoff, ls="--", alpha=.55)
ax.axhline(.1, color="gray", ls=":"); ax.set(xlabel="Pitch (nm)", ylabel="Contrast", ylim=(0, 1.05)); ax.legend()
fig.tight_layout(); fig.savefig(OUT/'resolution_curve.png',dpi=160); plt.close(fig)

pg = np.linspace(-1, 1, 301); pxx, pyy = np.meshgrid(pg, pg); pr = np.hypot(pxx, pyy); pt = np.arctan2(pyy, pxx)
pupil_cases = [("Defocus", ((2, 0, 1),)), ("Astigmatism", ((2, 2, 1),)),
               ("Coma", ((3, 1, 1),)), ("Spherical", ((4, 0, 1),))]
fig, axs = plt.subplots(1, 4, figsize=(12, 3))
for ax, (title, aberration) in zip(axs, pupil_cases):
    phase = np.where(pr <= 1, phase_map(pr, pt, aberration), np.nan)
    ax.imshow(phase, extent=(-1, 1, -1, 1), origin="lower", cmap="twilight")
    ax.set_title(title); ax.axis("off")
fig.tight_layout(); fig.savefig(OUT/'pupil_phases.png',dpi=160); plt.close(fig)
print(f"Wrote plots to {OUT}")
