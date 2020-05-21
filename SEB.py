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
df = df.resample("5min").mean()


dfT = pd.read_table("Desktop/IOP6_14oct/ES31Hz_IOP6_UTC.txt",sep = " ")
dfT["TIMESTAMP"]=pd.to_datetime(dfT["TIMESTAMP"])
dfT = dfT.set_index("TIMESTAMP") 
dfT = dfT.astype("float")
dfT = dfT.dropna()
dfT = dfT.resample("5min").mean()

dfS = pd.read_table("Desktop/IOP6_14oct/ES320Hz_IOP6_UTC.txt",sep = " ")
dfS["TIMESTAMP"]=pd.to_datetime(dfS["TIMESTAMP"])
dfS = dfS.set_index("TIMESTAMP") 
dfS = dfS.astype("float")
dfS = dfS.dropna()

dfF = pd.read_table("Desktop/IOP6_14oct/soilfluxes_IOP6_UTC.txt",sep = " ")
dfF["TIMESTAMP"]=pd.to_datetime(dfF["TIMESTAMP"])
dfF = dfF.set_index("TIMESTAMP") 
dfF = dfF.astype("float")
dfF = dfF.dropna()
dfF = dfF.resample("5min").mean()

#%% rotazione dati di vento e calcolo flussi turbolenti

uraw = dfS["u20_EM"]
vraw = dfS["v20_EM"]
wraw = dfS["w20_EM"]
Ts = dfS["ts20_EM"]
#[np,a]=size(uraw);

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
    
    um = np.nanmean(u)
    vm = np.nanmean(v)
    wm = np.nanmean(w)
    
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
    
    urm = np.nanmean(ur)
    vrm = np.nanmean(vr)
    wrm = np.nanmean(wr)
    
    nm = np.nanmean(bb+nx)
    
    urotm.append(urm)
    vrotm.append(vrm)
    wrotm.append(wrm)

    npm.append(nm)
    
dt = 36000

TsM = Ts.resample("30min").mean()

UP = []
VP = []
WP = []
TSP = []

WTP = []
UWP = []

WT = []
UW = []


for i in range(0,len(urotm)):
    uprimo = urot[i*dt:(i+1)*dt] - urotm[i]
    vprimo = vrot[i*dt:(i+1)*dt] - vrotm[i]
    wprimo = wrot[i*dt:(i+1)*dt] - wrotm[i]
    Tsprimo = Ts[i*dt:(i+1)*dt] - TsM[i]
    
    for j in range(0,len(uprimo)):
        WTP.append(uprimo[j]*Tsprimo[j])
        UWP.append(uprimo[j]*wprimo[j])


c = int(len(UWP)/6000.)
d = 6000

for i in range(0,c):
    WT.append(np.mean(WTP[i*d:(i+1)*d]))
    UW.append(np.mean(UWP[i*d:(i+1)*d]))



#%% R - net radiation 

LWinc = df["cg3up"]
LWrefl = df["cg3dn"]
SWinc = df["cm3up"]
SWrefl = df["cm3dn"]

R = np.zeros(len(LWinc))

for i in range(0,len(LWinc)):
    R[i] = SWinc[i] - SWrefl[i] + LWinc[i] - LWrefl[i]

#%% G - ground heat flux

G = dfF["Heat_flux_Avg"]

#%% H - sensible heat flux

md = 28.9647 # g/mol dry air molar mass
R = 8.314 # J/mol*K universal gas constant 
cp = 1000. # air specific heat  J/kg.K 
p = df["b_pressure"] # barometric pressure 
theta = dfT["Trh20_EM"]

H = np.zeros(len(theta))

for i in range(0,len(theta)-1):
    H[i] = (p[i]*0.01*md)/(R*(theta[i]+273.15))*cp*WT[i]
    
H[H<-500] = -300
#H[H>500] = 400
#    
#%% E - latent heat flux 
lambdaf = np.zeros(len(theta))

for i in range(0,len(theta)):
    lambdaf[i] = 3147.5-theta[i]*2.372
    
    
#%%
    
fig = plt.figure(figsize=(18, 5))
ax = fig.add_subplot(111)
ax.plot(np.arange(0,len(H)), H, label="H",marker = ".",linewidth=0.8,color="firebrick")
ax.plot(np.arange(0,len(G)), G, label="G",marker = ".",linewidth=0.8,color="tomato")
ax.plot(np.arange(0,len(R)), R, label="R",marker = ".",linewidth=0.8,color="royalblue")

plt.xlabel("time of the day (local time)",fontsize = "x-large")
plt.ylabel("SEB component",fontsize = "x-large")

plt.legend(loc='upper center', bbox_to_anchor=(1.08, 1.03));
ax.grid()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
plt.title("IOP6 - U (32m)",fontsize = "x-large")
plt.savefig("Desktop/SEB.png", quality = 95,dpi=300, bbox_inches = "tight")


    


















