"""
3-way comparison for TRACER Houston 2021-10-08 18 UTC:
  Fortran RRTM v3.3.1 (LW) and v2.7.2 (SW)
  Python buggy — rrtm_proj (original Python port)
  Python fixed — pyrrtm-lite (12 bugs corrected)
  ARM observation
"""
import sys, os, subprocess, shutil, tempfile, warnings
import numpy as np; import netCDF4 as nc
warnings.filterwarnings('ignore')

LW_BIN = '/home/vghate/Desktop/RRTM_LW/rrtm_v3.3.1_linux_pgf90'
SW_BIN = '/home/vghate/Desktop/rrtm_sw/rrtm_sw_linux_pgi_v2.7.2'
RULER1 = '0        1         2         3         4         5         6         7         8         9'
RULER2 = '123456789-' * 9

# ── Load TRACER sonde 18 UTC ──────────────────────────────────────────────────
f = nc.Dataset('/home/vghate/Desktop/TRACER_sonde/houinterpolatedsondeM1.c1.20211008.000030.nc')
z_km = np.array(f.variables['height'][:]); ti = 18*60
p_hPa = np.array(f.variables['bar_pres'][ti])*10.0
T_K   = np.array(f.variables['temp'][ti])+273.15
sh_gg = np.array(f.variables['sh'][ti]); f.close()

good = np.isfinite(p_hPa)&np.isfinite(T_K)&(p_hPa>0)
z_g=z_km[good][::-1]; p_g=p_hPa[good][::-1]   # surface first
T_g=T_K[good][::-1];  sh_g=sh_gg[good][::-1]
nlayers=len(z_g)-1

sys.path.insert(0, '/home/vghate/Desktop/rrtm_proj')
from rrtm_lw.std_atm import gas_vmr_std
wv_vmr = sh_g/(0.622+sh_g)
o3_vmr = gas_vmr_std(z_g[::-1],'o3')[::-1]

NA=6.02214076e23; G=9.80665; MAIR=28.97e-3
dp_Pa  = (p_g[:-1]-p_g[1:])*100.0
wbroad = dp_Pa*NA/(G*MAIR)/1e4
p_lay=0.5*(p_g[:-1]+p_g[1:]); T_lay=0.5*(T_g[:-1]+T_g[1:])
wv_lay=0.5*(wv_vmr[:-1]+wv_vmr[1:]); o3_lay=0.5*(o3_vmr[:-1]+o3_vmr[1:])
co2=422e-6; n2o=0.314e-6; ch4=1.9e-6; co=0.15e-6; o2=0.2095

print(f'Sonde: {nlayers} layers, P_sfc={p_g[0]:.1f} hPa, T_sfc={T_g[0]-273.15:.1f}°C')

def gas_row(wv, o3, wb):
    return f'{wv:15.7E}{co2:15.7E}{o3:15.7E}{n2o:15.7E}{co:15.7E}{ch4:15.7E}{o2:15.7E}{wb:15.7E}'

def lev_row(p_lav, t_lav, z_bot, p_bot, t_bot, z_top, p_top, t_top, first):
    s = f'{p_lav:15.7E}{t_lav:10.2f}{"":14}1'
    s += f'{z_bot:7.3f}{p_bot:8.2f}{t_bot:7.2f}' if first else ' '*22
    s += f'{z_top:7.3f}{p_top:8.3f}{t_top:7.2f}'
    return s

def write_lw(path):
    # Control line: match TROP example exactly
    # col 49='0'(IATM), col 69='0', col 83='0', col 87-90='  99'
    ctrl = ' '*49+'0'+' '*19+'0'+' '*13+'0'+'    99'
    with open(path,'w') as fh:
        fh.write('TRACER HOUSTON 2021-10-08 18UTC LW CLEAR-SKY CO2=422PPM\n')
        fh.write(RULER1+'\n'); fh.write(RULER2+'\n')
        fh.write('$ TRACER SOUNDING IATM=0 IOUT=99\n')
        fh.write(ctrl+'\n')
        fh.write(f'{T_g[0]:10.2f} 0\n')
        fh.write(f' 1{nlayers:3d}    7\n')
        for k in range(nlayers):
            fh.write(lev_row(p_lay[k],T_lay[k],
                             z_g[k],p_g[k],T_g[k],
                             z_g[k+1],p_g[k+1],T_g[k+1],k==0)+'\n')
            fh.write(gas_row(wv_lay[k],o3_lay[k],wbroad[k])+'\n')
        fh.write('%%%%%\n')

SZA=35.3; JULDAY=281

def write_sw(path):
    # SW control from working MLS example
    ctrl = ' '*49+'0'+' '*14+'0'+' 1  14   0  99'
    with open(path,'w') as fh:
        fh.write('TRACER HOUSTON 2021-10-08 18UTC SW CLEAR-SKY CO2=422PPM\n')
        fh.write(RULER1+'\n'); fh.write(RULER2+'\n')
        fh.write('$ TRACER SOUNDING IATM=0 IOSUN=0 IOUT=99\n')
        fh.write(ctrl+'\n')
        fh.write(f'              {JULDAY:6d}.000{SZA:10.3f}\n')
        fh.write('           2  '+''.join([f'{0.20:7.4f}' for _ in range(14)])+'\n')
        fh.write(f' 1{nlayers:3d}    7\n')
        for k in range(nlayers):
            fh.write(lev_row(p_lay[k],T_lay[k],
                             z_g[k],p_g[k],T_g[k],
                             z_g[k+1],p_g[k+1],T_g[k+1],k==0)+'\n')
            fh.write(gas_row(wv_lay[k],o3_lay[k],wbroad[k])+'\n')
        fh.write('%%%%%\n')

def run_rrtm(binary, inp, wdir):
    shutil.copy(inp, os.path.join(wdir,'INPUT_RRTM'))
    r = subprocess.run([binary], cwd=wdir, capture_output=True, text=True, timeout=60)
    out = os.path.join(wdir,'OUTPUT_RRTM')
    return open(out).read() if os.path.exists(out) else (r.stderr+r.stdout)

def parse_sfc(text, sw=False):
    rows=[]
    for l in text.splitlines():
        p=l.split()
        if len(p)>=(6 if sw else 4) and p[0].isdigit():
            try: rows.append((float(p[1]), float(p[5 if sw else 3]), float(p[2])))
            except: pass
    if not rows: return None,None,None
    rows.sort(key=lambda x:-x[0])
    return rows[0]

wd_lw=tempfile.mkdtemp(); wd_sw=tempfile.mkdtemp()
fi_lw=os.path.join(wd_lw,'i.txt'); fi_sw=os.path.join(wd_sw,'i.txt')
write_lw(fi_lw); write_sw(fi_sw)

print('Running Fortran RRTM_LW...')
lw_out=run_rrtm(LW_BIN,fi_lw,wd_lw)
print('Running Fortran RRTM_SW...')
sw_out=run_rrtm(SW_BIN,fi_sw,wd_sw)
shutil.rmtree(wd_lw); shutil.rmtree(wd_sw)

for label,out in [('LW',lw_out),('SW',sw_out)]:
    if 'FIO' in out or 'STOP' in out or 'error' in out.lower():
        print(f'{label} Fortran error:')
        for l in out.splitlines()[:6]: print(' ',l)

lw_r=parse_sfc(lw_out,sw=False); sw_r=parse_sfc(sw_out,sw=True)
fw_ld=lw_r[1]; fw_lu=lw_r[2]; fw_sd=sw_r[1]; fw_su=sw_r[2]

# ── Python runs ───────────────────────────────────────────────────────────────
sys.path.insert(0,'/home/vghate/Desktop/rrtm_sw_proj')
import rrtm_lw.run as lw_b_mod; from rrtm_lw.sounding import Sounding as SB
import rrtm_sw.run as sw_b_mod; from rrtm_sw.sounding import SWSounding as SWB
from rrtm_lw.std_atm import gas_vmr_std as gvmr_b

sys.path.insert(0,'/home/vghate/Desktop/pyrrtm-lite')
for k in list(sys.modules.keys()):
    if k in('lw','sw') or k.startswith('lw.') or k.startswith('sw.'): del sys.modules[k]
from lw.run import run as lw_f; from lw.sounding import Sounding as SF
from sw.run import run as sw_f; from sw.sounding import SWSounding as SWF
from lw.std_atm import gas_vmr_std as gvmr_f

zp=z_g[::-1]; pp=p_g[::-1]; Tp=T_g[::-1]-273.15; shp=sh_g[::-1]*1000
nlev=len(zp)

def gas(gvmr):
    return {'co2':np.full(nlev,422.),'o3':gvmr(zp,'o3')*1e6,
            'n2o':np.full(nlev,.314),'ch4':np.full(nlev,1.9),'co':np.full(nlev,.15)}

alb=np.full(14,0.20)
r_lb=lw_b_mod.run(sounding=SB(z_km=zp,T_C=Tp,P_hPa=pp,wv_gkg=shp,gas_ppm=gas(gvmr_b)),
                  coeffs_path='/home/vghate/Desktop/rrtm_proj/data/coeffs.nc')
r_lf=lw_f(sounding=SF(z_km=zp,T_C=Tp,P_hPa=pp,wv_gkg=shp,gas_ppm=gas(gvmr_f)),
           coeffs_path='/home/vghate/Desktop/pyrrtm-lite/data/lw_coeffs.nc')
r_sb=sw_b_mod.run(sounding=SWB(z_km=zp,T_C=Tp,P_hPa=pp,wv_gkg=shp,gas_ppm=gas(gvmr_b),sza_deg=SZA,albedo=alb),
                  coeffs_path='/home/vghate/Desktop/rrtm_sw_proj/data/coeffs.nc')
r_sf=sw_f(sounding=SWF(z_km=zp,T_C=Tp,P_hPa=pp,wv_gkg=shp,gas_ppm=gas(gvmr_f),sza_deg=SZA,albedo=alb),
          coeffs_path='/home/vghate/Desktop/pyrrtm-lite/data/sw_coeffs.nc')

lb_dn=float(r_lb['totdflux'][0]); lb_up=float(r_lb['totuflux'][-1])
lf_dn=float(r_lf['totdflux'][0]); lf_up=float(r_lf['totuflux'][-1])
sb_dn=float(r_sb['totdflux'][0]); sb_up=float(r_sb['totuflux'][-1])
sf_dn=float(r_sf['totdflux'][0]); sf_up=float(r_sf['totuflux'][-1])

# ARM obs 18 UTC
fobs=nc.Dataset('/home/vghate/Desktop/TRACER_sonde/houradflux1longM1.c2.20211008.060000.nc')
obs_lw=np.array(fobs.variables['downwelling_longwave'][:])
obs_sw=np.array(fobs.variables['downwelling_shortwave'][:])
bt=float(fobs.variables['base_time'][:]); toff=np.array(fobs.variables['time_offset'][:])
fobs.close()
ri=int(np.argmin(np.abs(bt+toff-(bt+18*3600-bt%86400))))
obs_lw18=float(obs_lw[ri]); obs_sw18=float(obs_sw[ri])

def fmt(v): return f'{v:9.2f}' if v is not None else '      n/a'

print()
print('='*78)
print('TRACER Houston 2021-10-08 18 UTC — Surface Downwelling Fluxes (W/m²)')
print(f'SZA={SZA}°  |  CO2=422 ppm  |  albedo=0.20  |  {nlayers} layers (0.008–25.5 km)')
print('='*78)
print(f'{"Code":32s}  {"LW↓ sfc":>9}  {"LW↑ TOA":>9}  {"SW↓ sfc":>9}  {"SW↑ TOA":>9}')
print('-'*78)
print(f'{"Fortran RRTM (v3.3.1 LW / v2.7.2 SW)":32s}  {fmt(fw_ld)}  {fmt(fw_lu)}  {fmt(fw_sd)}  {fmt(fw_su)}')
print(f'{"Python buggy — rrtm_proj":32s}  {lb_dn:9.2f}  {lb_up:9.2f}  {sb_dn:9.2f}  {sb_up:9.2f}')
print(f'{"Python fixed — pyrrtm-lite":32s}  {lf_dn:9.2f}  {lf_up:9.2f}  {sf_dn:9.2f}  {sf_up:9.2f}')
print(f'{"ARM observed (18 UTC)":32s}  {obs_lw18:9.2f}  {"n/a":>9}  {obs_sw18:9.2f}  {"n/a":>9}')
print('='*78)
print()
print('Differences vs Fortran RRTM:')
if fw_ld:
    print(f'  Py buggy − Fortran: LW↓={lb_dn-fw_ld:+.2f}  LW↑TOA={lb_up-fw_lu:+.2f}  SW↓={sb_dn-fw_sd:+.2f}  SW↑TOA={sb_up-fw_su:+.2f} W/m²')
    print(f'  Py fixed − Fortran: LW↓={lf_dn-fw_ld:+.2f}  LW↑TOA={lf_up-fw_lu:+.2f}  SW↓={sf_dn-fw_sd:+.2f}  SW↑TOA={sf_up-fw_su:+.2f} W/m²')
else:
    print('  Fortran output not available')
print(f'  Bug fix (fix−bug):  LW↓={lf_dn-lb_dn:+.2f}  LW↑TOA={lf_up-lb_up:+.2f}  SW↓={sf_dn-sb_dn:+.2f}  SW↑TOA={sf_up-sb_up:+.2f} W/m²')
print()
print('Differences vs ARM observation (18 UTC):')
if fw_ld:
    print(f'  Fortran RRTM − obs: LW↓={fw_ld-obs_lw18:+.2f}  SW↓={fw_sd-obs_sw18:+.2f} W/m²')
print(f'  Py buggy − obs:     LW↓={lb_dn-obs_lw18:+.2f}  SW↓={sb_dn-obs_sw18:+.2f} W/m²')
print(f'  Py fixed − obs:     LW↓={lf_dn-obs_lw18:+.2f}  SW↓={sf_dn-obs_sw18:+.2f} W/m²')
