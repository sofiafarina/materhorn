import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#%%

from datetime import datetime, timedelta

def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

              
# IOP6       
dts = [dt.strftime('%H:%M') for dt in 
       datetime_range(datetime(2012, 10, 14, 2), datetime(2012, 10, 15, 2), 
       timedelta(minutes=5))]
       
#%%
       
df = pd.read_table("Desktop/IOP6_14oct/IOP6_SEB_ES3_UTC.txt",sep = " ")
df["TIMESTAMP"]=pd.to_datetime(df["TIMESTAMP"])
df = df.set_index("TIMESTAMP") 
df = df.astype("float")
df = df.dropna()
#df = df.resample("5min").mean()


dfT = pd.read_table("Desktop/IOP6_14oct/ES31Hz_IOP6_UTC.txt",sep = " ")
dfT["TIMESTAMP"]=pd.to_datetime(dfT["TIMESTAMP"])
dfT = dfT.set_index("TIMESTAMP") 
dfT = dfT.astype("float")
dfT = dfT.dropna()
#dfT = dfT.resample("5min").mean()

dfS = pd.read_table("Desktop/IOP6_14oct/ES320Hz_IOP6_UTC.txt",sep = " ")
dfS["TIMESTAMP"]=pd.to_datetime(dfS["TIMESTAMP"])
dfS = dfS.set_index("TIMESTAMP") 
dfS = dfS.astype("float")
dfS = dfS.dropna()

dfSb = pd.read_table("Desktop/IOP6_14oct/ES320Hz_IOP6B_UTC.txt",sep = " ")
dfSb["TIMESTAMP"]=pd.to_datetime(dfSb["TIMESTAMP"])
dfSb = dfSb.set_index("TIMESTAMP") 
dfSb = dfSb.astype("float")
dfSb = dfSb.dropna()

dfF = pd.read_table("Desktop/IOP6_14oct/soilfluxes_IOP6_UTC.txt",sep = " ")
dfF["TIMESTAMP"]=pd.to_datetime(dfF["TIMESTAMP"])
dfF = dfF.set_index("TIMESTAMP") 
dfF = dfF.astype("float")
dfF = dfF.dropna()
dfF = dfF.resample("5min").mean()


#%% calcolo della w - mixing ratio

rhow = dfSb["rho_w"]
md = 28.9647 # g/mol dry air molar mass
Rc = 8.314 # J/mol*K universal gas constant 
p = df["b_pressure"] # barometric pressure 
theta = dfT["Trh20_EM"]+273.15

Wh = np.zeros(len(theta))

for i in range(0,len(theta)):
    Wh[i] = (Rc*(theta[i]))/(p[i]*md)*rhow[i]
    
    
#%% calcolo della w - mixing ratio

rhowM = dfSb["rho_w"].resample("30min").mean()
md = 28.9647 # g/mol dry air molar mass
Rc = 8.314 # J/mol*K universal gas constant 
pM = df["b_pressure"].resample("30min").mean() # barometric pressure 
thetaM = dfT["Trh20_EM"].resample("30min").mean()+273.15

WhM = np.zeros(len(thetaM))

for i in range(0,len(thetaM)):
    WhM[i] = (Rc*(thetaM[i]))/(pM[i]*md)*rhowM[i]
#%%
WhM = np.zeros(int(len(Wh)/1800))
for i in range(0,int(len(Wh)/1800)):
    WhM[i]=(np.mean(Wh[i*1800:(i+1)*1800]))


#%% rotazione dati di vento e calcolo flussi turbolenti

uraw = dfS["u20_EM"]
vraw = dfS["v20_EM"]
wraw = dfS["w20_EM"]
Ts = dfS["ts20_EM"]
rhow = dfSb["rho_w"]
rhowm = dfSb["rho_w"].resample("5min").mean()*0.001


n = len(uraw)
na_avg = 36000
urot = np.zeros(n)
vrot = np.zeros(n)
wrot = np.zeros(n)
urotm = []
vrotm = []
wrotm = []

npm = []


nx = np.arange(0,na_avg)

for bb in range(0,n-na_avg,na_avg):
    u = uraw.iloc[bb+nx]
    v = vraw.iloc[bb+nx]
    w = wraw.iloc[bb+nx]
    
    um = np.mean(u)
    vm = np.mean(v)
    wm = np.mean(w)
    
    s = np.sqrt(um**2+vm**2)
    theta = np.arctan2(vm,um)
    phi = np.arctan2(wm,s)
    
    rmatrix = np.zeros((3,3))
    rmatrix[0][0] = np.cos(phi)*np.cos(theta)
    rmatrix[0][1] = np.cos(phi)*np.sin(theta)
    rmatrix[0][2] = np.sin(phi)
    rmatrix[1][0] = -np.sin(theta)
    rmatrix[1][1] = np.cos(theta)
    rmatrix[1][2] = 0
    rmatrix[2][0] = -np.sin(phi)*np.cos(theta)
    rmatrix[2][1] = -np.sin(phi)*np.sin(theta)
    rmatrix[2][2] = np.cos(phi)
    
    
    X = np.vstack((u,v,w))
    XR = np.inner(rmatrix,X.T)
    ur = XR[0]
    vr = XR[1]
    wr = XR[2]
    
    for j in range(0,len(ur)-1):
        urot[bb+j] = ur[j]
        vrot[bb+j] = vr[j]
        wrot[bb+j] = wr[j]
    
    urm = np.mean(ur)
    vrm = np.mean(vr)
    wrm = np.mean(wr)
    
    nm = np.mean(bb+nx)
    
    urotm.append(urm)
    vrotm.append(vrm)
    wrotm.append(wrm)

    npm.append(nm)

# medie ogni cinque minuti 
u5min = np.zeros(288)
v5min = np.zeros(288)
w5min = np.zeros(288)

step = 6000

for i in range(0,288):
    u5min[i] = np.nanmean(urot[i*step:(i+1)*step])
    v5min[i] = np.nanmean(vrot[i*step:(i+1)*step])
    w5min[i] = np.nanmean(wrot[i*step:(i+1)*step])
    
# RECAP
# u5min media ogni cinque minuti
# urotm media ogni mezzora 
# urot 20Hz
    
dd = np.zeros(288)
for i in range(0,288):
    dd[i] = i*step
    
dt = 36000

TsM = Ts.resample("30min").mean()

UP = []
VP = []
WP = []
TSP = []
WhP = []

WTP = []
UWP = []
UUP = []
VVP = []
WWP = []
RWP = []
WWhP = []

WT = []
UW = []
UU = []
VV = []
WW = []
RW = []
WWh = []

WTmh = []
UWmh = []
UUmh = []
VVmh = []
WWmh = []
RWmh = []
WWhPmh = []


for i in range(0,len(urotm)):
    uprimo = np.zeros(dt)
    vprimo = np.zeros(dt)
    wprimo = np.zeros(dt)
    Tsprimo = np.zeros(dt)
    rhoprimo = np.zeros(dt)
    whprimo = np.zeros(dt)
    
    for k in range(0,dt):
        uprimo[k] = urot[i*dt+k] - urotm[i]
        vprimo[k] = vrot[i*dt+k] - vrotm[i]
        wprimo[k] = wrot[i*dt+k] - wrotm[i]
        Tsprimo[k] = Ts[i*dt+k] - TsM[i]
        rhoprimo[k] = rhow[i*dt+k] - rhowm[i]
        
        WTP.append(wprimo[k]*Tsprimo[k])
        UWP.append(uprimo[k]*wprimo[k])
        UUP.append(uprimo[k]*uprimo[k])
        VVP.append(vprimo[k]*vprimo[k])
        WWP.append(wprimo[k]*wprimo[k])
        RWP.append(wprimo[k]*rhoprimo[k])

c = int(len(UWP)/6000.)
d = 6000

for i in range(0,c):
    WT.append(np.mean(WTP[i*d:(i+1)*d]))
    UW.append(np.mean(UWP[i*d:(i+1)*d]))
    UU.append(np.mean(UUP[i*d:(i+1)*d]))
    VV.append(np.mean(VVP[i*d:(i+1)*d]))
    WW.append(np.mean(WWP[i*d:(i+1)*d]))
    RW.append(np.mean(RWP[i*d:(i+1)*d]))
    
    
#%% 
# calcolo delle fluttuazioni del mixing ratio del vapor d'acqua 
    
wrotave = np.zeros(len(df))
g = 20

for i in range(0,len(df)):
    wrotave[i] = np.mean(wrot[i*dt:(i+1)*dt])
    
ddt = 1800

for i in range(0,len(wrotm)-2):
    wprimo = np.zeros(ddt)
    whprimo = np.zeros(ddt)
    
    for k in range(0,ddt):
        wprimo[k] = wrotave[i*ddt+k] - wrotm[i]
        whprimo[k] = Wh[i*ddt+k] - WhM[i]
        WWhP.append(wprimo[k]*whprimo[k])


cc = int(len(WWhP)/300.)
dd = 300

for i in range(0,cc):
    WWh.append(np.mean(WWhP[i*dd:(i+1)*dd])) # ogni cinque min
    
cch = int(len(WWhP)/1800.)
ddmh = 1800

for i in range(0,cch):
    WWhPmh.append(np.mean(WTP[i*ddmh:(i+1)*ddmh])) # ogni mezzora

    
    
#%%
    

cmh = int(len(UWP)/36000.)
dmh = 36000

for i in range(0,cmh):
    WTmh.append(np.mean(WTP[i*dmh:(i+1)*dmh]))
    UWmh.append(np.mean(UWP[i*dmh:(i+1)*dmh]))
    UUmh.append(np.mean(UUP[i*dmh:(i+1)*dmh]))
    VVmh.append(np.mean(VVP[i*dmh:(i+1)*dmh]))
    WWmh.append(np.mean(WWP[i*dmh:(i+1)*dmh]))
    RWmh.append(np.mean(RWP[i*dmh:(i+1)*dmh]))
    
#%%
    
rad = df["NetTotR"].resample("5min").mean()


fig, ax1 = plt.subplots(figsize=(12, 5))

color = 'tab:blue'
ax1.set_xlabel('time')
ax1.set_ylabel("$\overline{u'u'}$", color=color)
ax1.plot(np.arange(0,len(UU)), UU)
ax1.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax1.set_ylim((0.0,5.5))
ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax2.set_ylabel('radiation [$W/m^2$]', color=color)  # we already handled the x-label with ax1
ax2.plot(np.arange(0,len(rad)), rad, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax2.xaxis.set_major_locator(plt.MaxNLocator(15))

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title("IOP 6, ES3, 32m: $\overline{u'u'}$")
plt.savefig("Desktop/uurad.png", quality = 95,dpi=300, bbox_inches = "tight")

#%%
    
rad = df["NetTotR"].resample("5min").mean()


fig, ax1 = plt.subplots(figsize=(12, 5))

color = 'tab:blue'
ax1.set_xlabel('time')
ax1.set_ylabel("$\overline{v'v'}$", color=color)
ax1.plot(np.arange(0,len(VV)), VV)
ax1.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax1.set_ylim((0.0,5.5))
ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax2.set_ylabel('radiation [$W/m^2$]', color=color)  # we already handled the x-label with ax1
ax2.plot(np.arange(0,len(rad)), rad, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax2.xaxis.set_major_locator(plt.MaxNLocator(15))

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title("IOP 6, ES3, 32m: $\overline{v'v'}$")
plt.savefig("Desktop/vvrad.png", quality = 95,dpi=300, bbox_inches = "tight")

#%%
    
rad = df["NetTotR"].resample("5min").mean()


fig, ax1 = plt.subplots(figsize=(12, 5))

color = 'tab:blue'
ax1.set_xlabel('time')
ax1.set_ylabel("$\overline{w'w'}$", color=color)
ax1.plot(np.arange(0,len(WW)), WW)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim((0.0,5.5))
fig.autofmt_xdate()
ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax2.set_ylabel('radiation [$W/m^2$]', color=color)  # we already handled the x-label with ax1
ax2.plot(np.arange(0,len(rad)), rad, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax2.xaxis.set_major_locator(plt.MaxNLocator(15))

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title("IOP 6, ES3, 32m: $\overline{w'w'}$")
plt.savefig("Desktop/wwrad.png", quality = 95,dpi=300, bbox_inches = "tight")


#%%
TKE = []

for i in range(0, len(UU)):
    TKE.append((UU[i]+VV[i]+WW[i])*0.5)

fig, ax1 = plt.subplots(figsize=(12, 5))

color = 'tab:blue'
ax1.set_xlabel('time')
ax1.set_ylabel("TKE", color=color)
ax1.plot(np.arange(0,len(TKE)), TKE)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim((0.0,5.5))
fig.autofmt_xdate()
ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax2.set_ylabel('radiation [$W/m^2$]', color=color)  # we already handled the x-label with ax1
ax2.plot(np.arange(0,len(rad)), rad, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax2.xaxis.set_major_locator(plt.MaxNLocator(15))

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title("IOP 6, ES3, 32m: TKE")
plt.savefig("Desktop/TKErad.png", quality = 95,dpi=300, bbox_inches = "tight")


#%%

fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(UU)), UU, label="Sonic Temp 20Hz",marker = ".",linewidth=0.8,color="firebrick")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("$\overline{u'u'}$",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - $\overline{u'u'}$",fontsize = "x-large")
plt.savefig("Desktop/uu.png", quality = 95,dpi=300, bbox_inches = "tight")

    
#%%

fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(Ts)), Ts, label="Sonic Temp 20Hz",marker = ".",linewidth=0.8,color="firebrick")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("Temperature (K)",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - Temperature 20Hz ",fontsize = "x-large")
plt.savefig("Desktop/TS.png", quality = 95,dpi=300, bbox_inches = "tight")

#%%

fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(TsM)), TsM, label="Sonic Temp 30min ave",marker = ".",linewidth=0.8,color="tomato")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("Temperature (K)",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - Temperature 30min ave",fontsize = "x-large")
plt.savefig("Desktop/TSave.png", quality = 95,dpi=300, bbox_inches = "tight")
    
#%%

fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(WTP)), WTP, label="$w'T'$ 20 Hz",marker = ".",linewidth=0.8,color="royalblue")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("$w'T'$",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - $w'T'$ ",fontsize = "x-large")
plt.savefig("Desktop/WTP.png", quality = 95,dpi=300, bbox_inches = "tight")
    
#%%

fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(WT)), WT, label="$w'T'$ 30 min ave",marker = ".",linewidth=0.8,color="tomato")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("$\overline{w'T'}$",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - $\overline{w'T'}$ ",fontsize = "x-large")
plt.savefig("Desktop/WT.png", quality = 95,dpi=300, bbox_inches = "tight")
    
#%%

fig, ax1 = plt.subplots(figsize=(12, 5))

color = 'tab:blue'
ax1.set_xlabel('time')
ax1.set_ylabel("$\overline{w'T'}$", color=color)
ax1.plot(np.arange(0,len(WT)), WT)
ax1.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
ax1.set_ylim((-3.5,3.0))
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax2.set_ylabel('radiation [$W/m^2$]', color=color)  # we already handled the x-label with ax1
ax2.plot(np.arange(0,len(rad)), rad, color=color)
ax2.tick_params(axis='y', labelcolor=color)
fig.autofmt_xdate()
ax2.xaxis.set_major_locator(plt.MaxNLocator(15))

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title("IOP 6, ES3, 32m: $\overline{w'T'}$")
plt.savefig("Desktop/wtrad.png", quality = 95,dpi=300, bbox_inches = "tight")


#%% R - net radiation 

LWinc = df["cg3up"]
LWrefl = df["cg3dn"]
SWinc = df["cm3up"]
SWrefl = df["cm3dn"]

R = np.zeros(len(LWinc))

for i in range(0,len(LWinc)):
    R[i] = SWinc[i] - SWrefl[i] + LWinc[i] - LWrefl[i]

#%% G - ground heat flux

G = -dfF["Heat_flux_Avg"]

#%% H - sensible heat flux

md = 28.9647 # g/mol dry air molar mass
Rc = 8.314 # J/mol*K universal gas constant 
cp = 1000. # air specific heat  J/kg.K 
p = df["b_pressure"].resample("5min").mean() # barometric pressure 
theta = dfT["Trh20_EM"].resample("5min").mean()+273.15

H = np.zeros(len(theta))

for i in range(0,len(theta)-1):
    H[i] = (p[i]*0.1*md)/(Rc*(theta[i]))*cp*WT[i]
    
cch = int(len(WT)/6.)
ddmh = 6
Hm = np.zeros(cch)
for i in range(0,cch):
    Hm[i] = (np.mean(H[i*ddmh:(i+1)*ddmh])) 
    
Hm[Hm<-500] = -300

Hm = - Hm
#H[H>500] = 400

#%% E - latent heat flux 

lambdaf = np.zeros(len(theta)-1)
E = np.zeros(len(theta))

for i in range(0,len(theta)-1):
    lambdaf[i] = 3147.5-theta[i]*2.372
    
for i in range(0,len(theta)-1):
    E[i] = lambdaf[i]*RW[i]*0.01
    
#%% E - latent heat flux  -- secondo tentativo
  
md = 28.9647 # g/mol dry air molar mass
Rc = 8.314 # J/mol*K universal gas constant 
cp = 1000. # air specific heat  J/kg.K 
p = df["b_pressure"].resample("30min").mean() # barometric pressure 
theta = dfT["Trh20_EM"].resample("30min").mean()+273.15

lambdaf = np.zeros(len(theta))
E = np.zeros(len(theta))

for i in range(0,len(theta)):
    lambdaf[i] = 3147.5-theta[i]*2.372

for i in range(0,len(theta)-3):
    E[i] = (p[i]*10*md)/(Rc*(theta[i]))*WWhPmh[i]
    
#%% TUTTO MEDIATO SU MEZZORA
    
LWinc = df["cg3up"].resample("30min").mean()
LWrefl = df["cg3dn"].resample("30min").mean()
SWinc = df["cm3up"].resample("30min").mean()
SWrefl = df["cm3dn"].resample("30min").mean()

R = np.zeros(len(LWinc))

for i in range(0,len(LWinc)):
    R[i] = SWinc[i] - SWrefl[i] + LWinc[i] - LWrefl[i]   

G = -dfF["Heat_flux_Avg"].resample("30min").mean()
    
#md = 28.9647 # g/mol dry air molar mass
#Rc = 8.314 # J/mol*K universal gas constant 
#cp = 1000. # air specific heat  J/kg.K 
#p = df["b_pressure"].resample("30min").mean() # barometric pressure 
#theta = dfT["Trh20_EM"].resample("30min").mean()+273.15
#
#
#
#H = np.zeros(len(theta))
#
#for i in range(0,len(theta)-1):
#    H[i] = (p[i]*0.1*md)/(Rc*(theta[i]))*cp*WT[i]
#    
#H[H<-500] = -400
    
md = 28.9647 # g/mol dry air molar mass
Rc = 8.314 # J/mol*K universal gas constant 
cp = 1000. # air specific heat  J/kg.K 
p = df["b_pressure"].resample("30min").mean() # barometric pressure 
theta = dfT["Trh20_EM"].resample("30min").mean()+273.15

lambdaf = np.zeros(len(theta))
E = np.zeros(len(theta))

for i in range(0,len(theta)):
    lambdaf[i] = 3147.5-theta[i]*2.372

for i in range(0,len(theta)-3):
    E[i] = (p[i]*100*md)/(Rc*(theta[i]))*WWhPmh[i]
    
    

    
fig = plt.figure(figsize=(15, 8))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(Hm)), Hm, label="H",marker = ".",linewidth=0.8,color="firebrick")
ax.plot(np.arange(0,len(G)), G, label="G",marker = ".",linewidth=0.8,color="tomato")
ax.plot(np.arange(0,len(R)), R, label="R",marker = ".",linewidth=0.8,color="royalblue")
ax.plot(np.arange(0,len(E)), E, label="E",marker = ".",linewidth=0.8,color="mediumturquoise")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("SEB component",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - Surface Energy Balance",fontsize = "x-large")
plt.savefig("Desktop/SEB.png", quality = 95,dpi=300, bbox_inches = "tight")


    

#%%
    
fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(WT)), WT, label="H",marker = ".",linewidth=0.8,color="firebrick")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("$\overline{w'T'}$",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - $\overline{w'T'}$",fontsize = "x-large")
plt.savefig("Desktop/WT.png", quality = 95,dpi=300, bbox_inches = "tight")


    
#%%
    
fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(UW)), UW, label="H",marker = ".",linewidth=0.8,color="royalblue")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("$\overline{w'u'}$",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - $\overline{w'u'}$",fontsize = "x-large")
plt.savefig("Desktop/UW.png", quality = 95,dpi=300, bbox_inches = "tight")


    



















