#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import matplotlib.pyplot as plt
import numpy as np
from euv_playground.multilayer import reflectivity_vs_angle, reflectivity_vs_wavelength
from euv_playground.zernike import zernike
from euv_playground.imaging import image_grating
from euv_playground.resist import plot_contrast_curve, shot_noise

OUT = ROOT/"plots"; OUT.mkdir(exist_ok=True)
lam=np.linspace(12.5,14.5,601); fig,ax=plt.subplots(); ax.plot(lam,reflectivity_vs_wavelength(lam)); ax.axvline(13.5,color='r',ls='--'); ax.set(xlabel='Wavelength (nm)',ylabel='Reflectivity',title='Mo/Si multilayer'); fig.tight_layout(); fig.savefig(OUT/'reflectivity_wavelength.png',dpi=160); plt.close(fig)
a=np.linspace(0,10,401); fig,ax=plt.subplots(); ax.plot(a,reflectivity_vs_angle(a)); ax.set(xlabel='Angle from normal (degrees)',ylabel='Reflectivity'); fig.tight_layout(); fig.savefig(OUT/'reflectivity_angle.png',dpi=160); plt.close(fig)
g=np.linspace(-1,1,301); xx,yy=np.meshgrid(g,g); rr=np.hypot(xx,yy); tt=np.arctan2(yy,xx); terms=[(2,0,'Defocus'),(2,2,'Astigmatism'),(2,-2,'Astigmatism 45°'),(3,1,'Coma X'),(3,-1,'Coma Y'),(3,3,'Trefoil'),(4,0,'Spherical'),(4,2,'Secondary astig.'),(4,4,'Quadrafoil')]; fig,axs=plt.subplots(3,3,figsize=(8,8));
for ax,(n,m,title) in zip(axs.flat,terms): ax.imshow(zernike(n,m,rr,tt),extent=(-1,1,-1,1),cmap='RdBu'); ax.set_title(title); ax.axis('off')
fig.tight_layout(); fig.savefig(OUT/'zernike_maps.png',dpi=160); plt.close(fig)
x=np.linspace(-1000,1000,4097); mask,clear=image_grating(x,64); _,blur=image_grating(x,64,zernike_terms=[(2,0,1.5)]); fig,ax=plt.subplots(); ax.plot(x,clear,label='No aberration'); ax.plot(x,blur,label='Defocus'); ax.set_xlim(-250,250); ax.set(xlabel='Position (nm)',ylabel='Intensity'); ax.legend(); fig.tight_layout(); fig.savefig(OUT/'imaging_demo.png',dpi=160); plt.close(fig)
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,4)); plot_contrast_curve(ax1); doses=np.linspace(5,100,200); _,p,n=shot_noise(doses,10); ax2.plot(doses,p,label='photons/pixel'); ax2.set(xlabel='Dose (mJ/cm²)',ylabel='Photons per 10 nm pixel'); axn=ax2.twinx(); axn.plot(doses,100*n,color='tab:red',label='noise'); axn.set_ylabel('Noise (%)',color='tab:red'); fig.tight_layout(); fig.savefig(OUT/'resist_and_shot_noise.png',dpi=160); plt.close(fig)
print(f"Wrote plots to {OUT}")
