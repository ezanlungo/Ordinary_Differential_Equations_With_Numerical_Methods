# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 13:13:48 2020

@author: EZEQUIEL
"""
import matplotlib.pyplot as plt
import numpy as np

#A da B
#B+C da A+C
#B+B da C

#A es x1 
#B es x2 
#C es x3

tmax=100.0
tmin=0.0
subst=np.array([10,100,500,1000])
h=(tmax-tmin)/subst
k1=0.04
k2=1e-3
k3=3e-3
x1=1
x2=0
x3=0


def xprima (x1,x2,x3):
    A=-k1*x1+k2*x2*x3
    B=k1*x1-k2*x2*x3-k3*x2**2
    C=k3*x2**2
    return A,B,C

euler_p_error=np.array([])
eulera=[]
eulerb=[]
eulerc=[]
for j in range (len(h)):
    t=np.linspace(tmin,tmax,subst[j]+1)
    euler_a=np.array([])
    euler_b=np.array([])
    euler_c=np.array([])
    euler_a=np.append(euler_a,x1)
    euler_b=np.append(euler_b,x2)
    euler_c=np.append(euler_c,x3)
    for i in range (len(t)-1):
        a=euler_a[i]+h[j]*xprima(euler_a[i],euler_b[i],euler_c[i])[0]
        b=euler_b[i]+h[j]*xprima(euler_a[i],euler_b[i],euler_c[i])[1]
        c=euler_c[i]+h[j]*xprima(euler_a[i],euler_b[i],euler_c[i])[2]
        euler_a=np.append(euler_a,a)
        euler_b=np.append(euler_b,b)
        euler_c=np.append(euler_c,c)
    eulera.append(euler_a)
    eulerb.append(euler_b)
    eulerc.append(euler_c)
    euler_p_error=np.append(euler_p_error,euler_b[-1])


taylor_p_error=np.array([])
taylora=[]
taylorb=[]
taylorc=[]
for j in range (len(h)):
    t=np.linspace(tmin,tmax,subst[j]+1)
    taylor_a=np.array([])
    taylor_b=np.array([])
    taylor_c=np.array([])    
    taylor_a=np.append(taylor_a,x1)
    taylor_b=np.append(taylor_b,x2)
    taylor_c=np.append(taylor_c,x3)    
    for i in range (len(t)-1):
        xsegunda=-k1*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[0]+k2*taylor_c[i]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[1]+k2*taylor_b[i]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[2]
        ysegunda=k1*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[0]-k2*taylor_c[i]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[1]+k2*taylor_b[i]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[2]-2*k3*taylor_b[i]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[1]
        zsegunda=2*k3*taylor_b[i]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[1] 
        a=taylor_a[i]+h[j]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[0]+((h[j]**2)/2)*xsegunda
        b=taylor_b[i]+h[j]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[1]+((h[j]**2)/2)*ysegunda
        c=taylor_c[i]+h[j]*xprima(taylor_a[i],taylor_b[i],taylor_c[i])[2]+((h[j]**2)/2)*zsegunda
        taylor_a=np.append(taylor_a,a)
        taylor_b=np.append(taylor_b,b)
        taylor_c=np.append(taylor_c,c)
    taylora.append(taylor_a)
    taylorb.append(taylor_b)
    taylorc.append(taylor_c)
    taylor_p_error=np.append(taylor_p_error,taylor_b[-1])

 
rk_p_error=np.array([])
rka=[]
rkb=[]
rkc=[]
for j in range (len(h)):
    t=np.linspace(tmin,tmax,subst[j]+1)
    rk_a=np.array([])
    rk_b=np.array([])
    rk_c=np.array([])
    rk_a=np.append(rk_a,x1)
    rk_b=np.append(rk_b,x2)
    rk_c=np.append(rk_c,x3)
    for i in range(len(t)-1):
        k1a=xprima(rk_a[i],rk_b[i],rk_c[i])[0]
        k1b=xprima(rk_a[i],rk_b[i],rk_c[i])[1]
        k1c=xprima(rk_a[i],rk_b[i],rk_c[i])[2]
        k2a=xprima(rk_a[i]+1/2.0*h[j]*k1a,rk_b[i]+1/2.0*h[j]*k1b,rk_c[i]+1/2.0*h[j]*k1c)[0]
        k2b=xprima(rk_a[i]+1/2.0*h[j]*k1a,rk_b[i]+1/2.0*h[j]*k1b,rk_c[i]+1/2.0*h[j]*k1c)[1]
        k2c=xprima(rk_a[i]+1/2.0*h[j]*k1a,rk_b[i]+1/2.0*h[j]*k1b,rk_c[i]+1/2.0*h[j]*k1c)[2]
        k3a=xprima(rk_a[i]+1/2.0*h[j]*k2a,rk_b[i]+1/2.0*h[j]*k2b,rk_c[i]+1/2.0*h[j]*k2c)[0]
        k3b=xprima(rk_a[i]+1/2.0*h[j]*k2a,rk_b[i]+1/2.0*h[j]*k2b,rk_c[i]+1/2.0*h[j]*k2c)[1]
        k3c=xprima(rk_a[i]+1/2.0*h[j]*k2a,rk_b[i]+1/2.0*h[j]*k2b,rk_c[i]+1/2.0*h[j]*k2c)[2]  
        k4a=xprima(rk_a[i]+h[j]*k3a,rk_b[i]+h[j]*k3b,rk_c[i]+h[j]*k3c)[0]
        k4b=xprima(rk_a[i]+h[j]*k3a,rk_b[i]+h[j]*k3b,rk_c[i]+h[j]*k3c)[1]
        k4c=xprima(rk_a[i]+h[j]*k3a,rk_b[i]+h[j]*k3b,rk_c[i]+h[j]*k3c)[2]
        rk_a=np.append(rk_a,rk_a[i]+h[j]/6.0*(k1a+2*k2a+3*k3a+k4a))       
        rk_b=np.append(rk_b,rk_b[i]+h[j]/6.0*(k1b+2*k2b+3*k3b+k4b)) 
        rk_c=np.append(rk_c,rk_c[i]+h[j]/6.0*(k1c+2*k2c+3*k3c+k4c))
    rka.append(rk_a)
    rkb.append(rk_b)
    rkc.append(rk_c)
    rk_p_error=np.append(rk_p_error,rk_b[-1])


ab_p_error=np.array([])
aba=[]
abb=[]
abc=[] 
for j in range (len(h)):
    t=np.linspace(tmin,tmax,subst[j]+1)
    ab_a=np.array([])
    ab_b=np.array([])
    ab_c=np.array([])    
    ab_a=np.append(ab_a,x1)
    ab_b=np.append(ab_b,x2)
    ab_c=np.append(ab_c,x3) 
    for i in range(len(t)-1): 
        if i<3:
            k1a2=xprima(ab_a[i],ab_b[i],ab_c[i])[0]
            k1b2=xprima(ab_a[i],ab_b[i],ab_c[i])[1]
            k1c2=xprima(ab_a[i],ab_b[i],ab_c[i])[2]
            k2a2=xprima(ab_a[i]+1/2.0*h[j]*k1a2,ab_b[i]+1/2.0*h[j]*k1b2,ab_c[i]+1/2.0*h[j]*k1c2)[0]
            k2b2=xprima(ab_a[i]+1/2.0*h[j]*k1a2,ab_b[i]+1/2.0*h[j]*k1b2,ab_c[i]+1/2.0*h[j]*k1c2)[1]
            k2c2=xprima(ab_a[i]+1/2.0*h[j]*k1a2,ab_b[i]+1/2.0*h[j]*k1b2,ab_c[i]+1/2.0*h[j]*k1c2)[2]
            k3a2=xprima(ab_a[i]+1/2.0*h[j]*k2a2,ab_b[i]+1/2.0*h[j]*k2b2,ab_c[i]+1/2.0*h[j]*k2c2)[0]
            k3b2=xprima(ab_a[i]+1/2.0*h[j]*k2a2,ab_b[i]+1/2.0*h[j]*k2b2,ab_c[i]+1/2.0*h[j]*k2c2)[1]
            k3c2=xprima(ab_a[i]+1/2.0*h[j]*k2a2,ab_b[i]+1/2.0*h[j]*k2b2,ab_c[i]+1/2.0*h[j]*k2c2)[2]  
            k4a2=xprima(ab_a[i]+h[j]*k3a2,rk_b[i]+h[j]*k3b2,rk_c[i]+h[j]*k3c2)[0]
            k4b2=xprima(ab_a[i]+h[j]*k3a2,rk_b[i]+h[j]*k3b2,rk_c[i]+h[j]*k3c2)[1]
            k4c2=xprima(ab_a[i]+h[j]*k3a2,rk_b[i]+h[j]*k3b2,rk_c[i]+h[j]*k3c2)[2]
            ab_a=np.append(ab_a,ab_a[i]+h[j]/6.0*(k1a2+2*k2a2+3*k3a2+k4a2))       
            ab_b=np.append(ab_b,ab_b[i]+h[j]/6.0*(k1b2+2*k2b2+3*k3b2+k4b2)) 
            ab_c=np.append(ab_c,ab_c[i]+h[j]/6.0*(k1c2+2*k2c2+3*k3c2+k4c2))
        else:
            pred0a=xprima(ab_a[i],ab_b[i],ab_c[i])[0]
            pred0b=xprima(ab_a[i],ab_b[i],ab_c[i])[1]
            pred0c=xprima(ab_a[i],ab_b[i],ab_c[i])[2]
            pred1a=xprima(ab_a[i-1],ab_b[i-1],ab_c[i-1])[0]
            pred1b=xprima(ab_a[i-1],ab_b[i-1],ab_c[i-1])[1]
            pred1c=xprima(ab_a[i-1],ab_b[i-1],ab_c[i-1])[2]
            pred2a=xprima(ab_a[i-2],ab_b[i-2],ab_c[i-2])[0]
            pred2b=xprima(ab_a[i-2],ab_b[i-2],ab_c[i-2])[1]
            pred2c=xprima(ab_a[i-2],ab_b[i-2],ab_c[i-2])[2]
            pred3a=xprima(ab_a[i-3],ab_b[i-3],ab_c[i-3])[0]
            pred3b=xprima(ab_a[i-3],ab_b[i-3],ab_c[i-3])[1]
            pred3c=xprima(ab_a[i-3],ab_b[i-3],ab_c[i-3])[2]
            preda=ab_a[i]+h[j]/24.0*(55*pred0a-59*pred1a+37*pred2a-9*pred3a)  
            predb=ab_b[i]+h[j]/24.0*(55*pred0b-59*pred1b+37*pred2b-9*pred3b) 
            predc=ab_c[i]+h[j]/24.0*(55*pred0c-59*pred1c+37*pred2c-9*pred3c) 
            
            corr1a=xprima(preda,predb,predc)[0]
            corr1b=xprima(preda,predb,predc)[1]
            corr1c=xprima(preda,predb,predc)[2]
            corra=ab_a[i]+h[j]/24.0*(9*corr1a+19*pred0a-5*pred1a+pred2a)
            corrb=ab_b[i]+h[j]/24.0*(9*corr1b+19*pred0b-5*pred1b+pred2b)
            corrc=ab_c[i]+h[j]/24.0*(9*corr1c+19*pred0c-5*pred1c+pred2c)
            ab_a=np.append(ab_a,corra)
            ab_b=np.append(ab_b,corrb)
            ab_c=np.append(ab_c,corrc)
    aba.append(ab_a)
    abb.append(ab_b)
    abc.append(ab_c)
    ab_p_error=np.append(ab_p_error,ab_b[-1])


substanalitico=20000
hanalitico=(tmax-tmin)/substanalitico
tanalitico=np.linspace(tmin,tmax,substanalitico+1)
    
analitico_a=np.array([])
analitico_b=np.array([])
analitico_c=np.array([])

analitico_a=np.append(analitico_a,x1)
analitico_b=np.append(analitico_b,x2)
analitico_c=np.append(analitico_c,x3) 

for i in range(len(tanalitico)-1):
    k1aa=xprima(analitico_a[i],analitico_b[i],analitico_c[i])[0]
    k1bb=xprima(analitico_a[i],analitico_b[i],analitico_c[i])[1]
    k1cc=xprima(analitico_a[i],analitico_b[i],analitico_c[i])[2]
    k2aa=xprima(analitico_a[i]+1/2.0*hanalitico*k1aa,analitico_b[i]+1/2.0*hanalitico*k1bb,analitico_c[i]+1/2.0*hanalitico*k1cc)[0]
    k2bb=xprima(analitico_a[i]+1/2.0*hanalitico*k1aa,analitico_b[i]+1/2.0*hanalitico*k1bb,analitico_c[i]+1/2.0*hanalitico*k1cc)[1]
    k2cc=xprima(analitico_a[i]+1/2.0*hanalitico*k1aa,analitico_b[i]+1/2.0*hanalitico*k1bb,analitico_c[i]+1/2.0*hanalitico*k1cc)[2]
    k3aa=xprima(analitico_a[i]+1/2.0*hanalitico*k2aa,analitico_b[i]+1/2.0*hanalitico*k2bb,analitico_c[i]+1/2.0*hanalitico*k2cc)[0]
    k3bb=xprima(analitico_a[i]+1/2.0*hanalitico*k2aa,analitico_b[i]+1/2.0*hanalitico*k2bb,analitico_c[i]+1/2.0*hanalitico*k2cc)[1]
    k3cc=xprima(analitico_a[i]+1/2.0*hanalitico*k2aa,analitico_b[i]+1/2.0*hanalitico*k2bb,analitico_c[i]+1/2.0*hanalitico*k2cc)[2]  
    k4aa=xprima(analitico_a[i]+hanalitico*k3aa,analitico_b[i]+hanalitico*k3bb,analitico_c[i]+hanalitico*k3cc)[0]
    k4bb=xprima(analitico_a[i]+hanalitico*k3aa,analitico_b[i]+hanalitico*k3bb,analitico_c[i]+hanalitico*k3cc)[1]
    k4cc=xprima(analitico_a[i]+hanalitico*k3aa,analitico_b[i]+hanalitico*k3bb,analitico_c[i]+hanalitico*k3cc)[2]
    analitico_a=np.append(analitico_a,analitico_a[i]+hanalitico/6.0*(k1aa+2*k2aa+3*k3aa+k4aa))       
    analitico_b=np.append(analitico_b,analitico_b[i]+hanalitico/6.0*(k1bb+2*k2bb+3*k3bb+k4bb)) 
    analitico_c=np.append(analitico_c,analitico_c[i]+hanalitico/6.0*(k1cc+2*k2cc+3*k3cc+k4cc))


fig, axs = plt.subplots(nrows=1,ncols=len(h),sharey=True)
for p in range (len(h)):
    t=np.linspace(tmin,tmax,subst[p]+1)
    axs[p].set_title('Euler con %i subdivisiones' %subst[p])
    axs[p].plot(t,eulera[p],'r',label='Ca')
    axs[p].plot(t,eulerb[p],'b',label='Cb')
    axs[p].plot(t,eulerc[p],'g',label='Cc')
    axs[p].set_xlabel("Tiempo")
    axs[0].set_ylabel("C [mol/L]")
    axs[p].grid(True)
    axs[p].legend(loc='best')


fig, axs = plt.subplots(nrows=1,ncols=len(h),sharey=True)
for p in range (len(h)):
    t=np.linspace(tmin,tmax,subst[p]+1)
    axs[p].set_title('Taylor con %i subdivisiones' %subst[p])
    axs[p].plot(t,taylora[p],'r',label='Ca')
    axs[p].plot(t,taylorb[p],'b',label='Cb')
    axs[p].plot(t,taylorc[p],'g',label='Cc')
    axs[p].set_xlabel("Tiempo")
    axs[0].set_ylabel("C [mol/L]")
    axs[p].grid(True)
    axs[p].legend(loc='best')
    

fig, axs = plt.subplots(nrows=1,ncols=len(h),sharey=True)
for p in range (len(h)):
    t=np.linspace(tmin,tmax,subst[p]+1)
    axs[p].set_title('RK 4 con %i subdivisiones' %subst[p])
    axs[p].plot(t,rka[p],'r',label='Ca')
    axs[p].plot(t,rkb[p],'b',label='Cb')
    axs[p].plot(t,rkc[p],'g',label='Cc')
    axs[p].set_xlabel("Tiempo")
    axs[0].set_ylabel("C [mol/L]")
    axs[p].grid(True)
    axs[p].legend(loc='best')


fig, axs = plt.subplots(nrows=1,ncols=len(h),sharey=True)
for p in range (len(h)):
    t=np.linspace(tmin,tmax,subst[p]+1)
    axs[p].set_title('P-C con %i subdivisiones' %subst[p])
    axs[p].plot(t,aba[p],'r',label='Ca')
    axs[p].plot(t,abb[p],'b',label='Cb')
    axs[p].plot(t,abc[p],'g',label='Cc')
    axs[p].set_xlabel("Tiempo")
    axs[0].set_ylabel("C [mol/L]")
    axs[p].grid(True)
    axs[p].legend(loc='best')    

fig5 = plt.figure()
plt.title('Comparación')
plt.plot(t,eulerb[-1],'r--',label='Cb Euler')    
plt.plot(t,taylorb[-1],'b--',label='Cb Taylor')    
plt.plot(t,rkb[-1],'g--',label='Cb RK4')    
plt.plot(t,abb[-1],'y--',label='Cb P-C')    
plt.plot(tanalitico,analitico_b,'k-',label='Cb Solución')    
plt.xlabel("Tiempo") 
plt.ylabel("Concentración [mol/L]")
plt.grid(True) 
plt.legend(loc='center right')
plt.show()

fig6 = plt.figure()
plt.title('Errores según tamaño de paso para Cb')
plt.plot(h,np.abs(analitico_b[-1]-euler_p_error), 'r-',label='Euler')
plt.plot(h,np.abs(analitico_b[-1]-taylor_p_error), 'b-',label='Taylor')
plt.plot(h,np.abs(analitico_b[-1]-rk_p_error), 'g-',label='RK4')
plt.plot(h,np.abs(analitico_b[-1]-ab_p_error), 'y-',label='Predictor-Corrector')
plt.xlabel("Tamaño de paso") 
plt.ylabel("Error Global")
plt.grid(True)
plt.legend(loc='center right')
plt.show()

for q in range (len(h)):
    print ('\nError para el cálculo de la concentración de B con',np.round(subst[q]),'subdivisiones: \n')
    print ('El error global para Euler es:', np.abs(np.round(euler_p_error[q]-analitico_b[-1],10)),)
    print ('El error global para Taylor es:', np.abs(np.round(taylor_p_error[q]-analitico_b[-1],10)),)
    print ('El error global para Runge-Kutta 4 es:', np.abs(np.round(rk_p_error[q]-analitico_b[-1],10)),)
    print ('El error global para el Método P-C es:', np.abs(np.round(ab_p_error[q]-analitico_b[-1],10)),)