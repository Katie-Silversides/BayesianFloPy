# -*- coding: utf-8 -*-
"""
Uses Bayesian Model selection / MCMC to calculate the variables that contribute to 
aquifer inputs and outputs.
Comparison using head values of monitoring wells
"""
import numpy as np
from itertools import product
import os
import time
import flopy



np.random.seed(1)
# os.chdir('.')
out_folder='./output_q_900-1100_r_005-01'
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

var_cons=np.loadtxt('./parameters_FP_ModelS_1.csv', delimiter=',', skiprows=1, usecols=[1,2])

nu=2 # number of uncertain components
model_nu = 2**nu #number of model options
theta=np.array(list(product([0,1],repeat=nu)))

blocksize=1000 # runs with 1000 parameter values at the time
runs=2
fail_n = 0 #number of times modflow fails to terminate normally
complete_n = 0 #number of times modflow terminates normally

sd=0.025 # standard deviation for relative error
m=np.arange(0,len(theta))

# import correct monitoring wells
t_well = np.loadtxt('FloPy_Test1b_monitoring_wells.csv',delimiter=',')

#%%

#%%
start = time.time()

li = np.empty([model_nu,blocksize])

for run in range(0,runs):
    # sample random points between 0 and 1
    a=np.random.uniform(size=[len(var_cons),blocksize]) #variables

    for i in range(0,len(var_cons)):
        #scale random points to parameter constraints in var_cons
        # (maximum-minimum)*x+minimum
        a[i,:]=(var_cons[i,1]-var_cons[i,0])*a[i,:]+var_cons[i,0]
    
    for block in range(blocksize): # for each block of models in the run
    
        b_vars = a[:,block]
        
        for var in range(model_nu): #for each model option in the block
            
        # Get the variables for this model
            m_vars = b_vars * theta[var,:]
            
            # Set up FloPy model - these are the parts that don't change

            ## Set variables
            # specified head
            h1 = 100 

            #Number of layers
            Nlay = 10

            #Number of rows and columns
            N = 101

            #lengths of the sides of the model
            L = 400.0

            #Aquifer thickness
            H = 50.0

            #Hydrolic conductivity
            k = 1.0

            mf6_path = "C:\\Users\\ksil8584\\OneDrive - The University of Sydney (Staff)\\Documents\\Modflow\mf6.4.1\\bin\\mf6.exe"
            name = "FloPy_ModelS_T1b"
            sim = flopy.mf6.MFSimulation(
                sim_name=name, exe_name=mf6_path, version="mf6", sim_ws=out_folder
            )

            ## Create FloPy temporal discretization (TDIS) object
            tdis = flopy.mf6.ModflowTdis(
                sim, pname="tdis", time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)] # add  for _ in range(2)  to persiod data if nper >1
            )

            ## Create FloPy iterative model solution (IMS) object
            ims = flopy.mf6.ModflowIms(
                sim,
                pname="ims",
                complexity="SIMPLE",
                linear_acceleration="BICGSTAB",
            )

            ## Create Flopy groundwater flow (gwf) model object
            model_nam_file = f"{name}.nam"
            gwf = flopy.mf6.ModflowGwf(
                sim,
                modelname=name,
                model_nam_file=model_nam_file,
                save_flows=True,
                newtonoptions="NEWTON UNDER_RELAXATION",
            )

            ## Create discretization (DIS) Package
            # All layers are given equal thickness
            bot = np.linspace(-H / Nlay, -H, Nlay)
            #  delrow and delcol are computed from model size L and number of cells
            delrow = delcol = L / (N - 1)
            dis = flopy.mf6.ModflowGwfdis(
                gwf,
                nlay=Nlay,
                nrow=N,
                ncol=N,
                delr=delrow,
                delc=delcol,
                top=0.0,
                botm=bot,
            )

            ## Create initial conditions (IC) Package
            start = h1 * np.ones((Nlay, N, N))
            ic = flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=start)

            ## Create node property flow (NPF) package
            npf = flopy.mf6.ModflowGwfnpf(
                gwf,
                icelltype=1,
                k=k,
            )
            
            ## Create storage (STO) package
            sto = flopy.mf6.ModflowGwfsto(
                gwf, 
                iconvert = 1
            )

            ## Create constant head (CHD) Package
            chd_rec = []
            layer = 0
            for row_col in range(0, N):
                chd_rec.append(((layer, row_col, 0), h1))
                chd_rec.append(((layer, row_col, N - 1), h1))
                if row_col != 0 and row_col != N - 1:
                    chd_rec.append(((layer, 0, row_col), h1))
                    chd_rec.append(((layer, N - 1, row_col), h1))
            chd = flopy.mf6.ModflowGwfchd(
                gwf,
                stress_period_data=chd_rec,
            )

            ## Create output control (OC) Package
            headfile = f"{name}.hds"
            head_filerecord = [headfile]
            budgetfile = f"{name}.cbb"
            budget_filerecord = [budgetfile]
            saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
            printrecord = [("HEAD", "LAST")]
            oc = flopy.mf6.ModflowGwfoc(
                gwf,
                saverecord=saverecord,
                head_filerecord=head_filerecord,
                budget_filerecord=budget_filerecord,
                printrecord=printrecord,
            )

            
            # Run Modflow model using variables that change
            #Well pumping rate
            q = -m_vars[0]
    
            ## Create well (WEL) package
            wel_rec = [(Nlay - 1, int(N / 4), int(N / 4), q)]
            wel = flopy.mf6.ModflowGwfwel(
                gwf,
                stress_period_data=wel_rec,
            )
    
            # create recharge package
            # recharge rate
            rch_rate=[m_vars[1]]
            
            #Creach recharge (rcha) package
            rch = flopy.mf6.ModflowGwfrcha(
                gwf, 
                readasarrays=True, 
                pname="rch", 
                recharge=rch_rate
            )
            
            ## Write the datasets
            sim.write_simulation()

            ## Run simulation
            try:
                sim_tru = sim.run_simulation()
            except:
                print("MODFLOW 6 did not terminate normally.")
                
            if sim_tru[0]==True:
                # create monitoring wells
                h = gwf.output.head().get_data(kstpkper=(0, 0))
                m_well = np.empty([4,4])
                for m_x in range(20,100,20):
                    m_xi = int(m_x/20-1)
                    for m_y in range(20,100,20):
                        m_yi = int(m_y/20-1)
                        m_well[m_xi,m_yi] = h[0][m_x,m_y]
                complete_n += 1
            else:
                m_well = np.full([4,4],150)
                fail_n += 1
            
                    
            #compare monitoring wells to true values
            well_dif = abs(t_well-m_well)
    
            # Calculate average height error            
            delta=np.mean(well_dif)
            #calculate average height
            well_ht = np.mean(t_well)
            # Calculate relative error
            delta_r=delta/well_ht
            # Calculate likelihood based on relative water balance error
            li[var,block]=1/sd*np.exp((-(1/2*delta_r**2)/sd**2))
            
            np.save(os.path.join(out_folder, 'well_dif'+str(run)+str(block)+'.npy'),well_dif)
            np.save(os.path.join(out_folder, 'm_well'+str(run)+str(block)+'.npy'),m_well)


# Metropolis-Hastings sampling on likelihood matrix
    accept=np.zeros(shape=[blocksize,1],dtype=bool)
    models=np.zeros(shape=[blocksize,1],dtype=int)
    apost=np.zeros(shape=[2,blocksize], dtype=float)
    for p in np.arange(0,blocksize): #stepping through likelihood matrix
        q=np.random.choice(m, p=li[:,p]/np.sum(li[:,p])) #choose model based on likelihood
        prop=li[q,p] # proposal likelihood
        if p==0 and run==0: #Only starts from zero if.  
            current=prop # current likelihood
        else:
            accept[p]=np.random.rand() <= prop / current
            
            if accept[p]==True: # accept new state
                current=prop  # current likelihood updated
                models[p]=q
                apost[:,p]=a[:,p]
            elif accept[p]==False: # reject and move old state forward
                apost[:,p]=apost[:,p-1]
                models[p]=models[p-1]
                    
    np.save(os.path.join(out_folder, 'accept'+str(run)+'.npy'),accept)
    np.save(os.path.join(out_folder, 'models'+str(run)+'.npy'),models)
    np.save(os.path.join(out_folder, 'vars_po'+str(run)+'.npy'),apost)
    np.save(os.path.join(out_folder, 'vars_pr'+str(run)+'.npy'),a)
    np.save(os.path.join(out_folder, 'delta_r'+str(run)+'.npy'),delta_r)
    np.save(os.path.join(out_folder, 'li'+str(run)+'.npy'),li)
    
np.save(os.path.join(out_folder, 'runs.npy'),runs)
np.save(os.path.join(out_folder, 'blocksize.npy'),blocksize)
np.save(os.path.join(out_folder, 'variable_constraints.npy'), var_cons)

#print('It took {0:0.1f} seconds'.format(time.time() - start))
del(apost)
del(a)
del(delta_r)
del(delta)
#%%
accept=np.zeros(shape=[runs*blocksize,1],dtype=bool)

for run in np.arange(0,runs):    
    accept[run*blocksize:(1+run)*blocksize]=np.load(os.path.join(out_folder, 
           'accept'+str(run)+'.npy'))
accept=np.squeeze(accept)    

models = np.zeros(shape=[runs*blocksize,1],dtype=int)   
for run in range(0,runs):    
    models[run*blocksize:(1+run)*blocksize]=np.load(os.path.join(out_folder, 
           'models'+str(run)+'.npy'))
models=np.squeeze(models)
np.save(os.path.join(out_folder, 'models_vector.npy'),models)
model_sum=np.zeros(shape=[2**nu,1],dtype=int)    

for num in range(0,2**nu):
    model_sum[num]=np.sum(models==num)
np.save(os.path.join(out_folder, 'model_sum.npy'),model_sum)

# calc results
sim_res = np.sum(model_sum*theta,axis=0)
np.save(os.path.join(out_folder, 'sim_res.npy'),sim_res)

sim_res_p = sim_res/sum(model_sum)
np.save(os.path.join(out_folder, 'sim_res_p.npy'),sim_res_p)

cond_res_t = np.empty([nu,nu])
cond_res_f = np.empty([nu,nu])
cond_res_t_p = np.empty([nu,nu])
cond_res_f_p = np.empty([nu,nu])
model_theta = model_sum*theta

for idx in range(nu):
    model_t = model_theta[theta[:,idx]==1]
    model_sum_t = model_sum[theta[:,idx]==1]
    cond_res_t[:,idx] = np.sum(model_t,axis=0)
    cond_res_t_p[:,idx] = (np.sum(model_t,axis=0))/np.sum(model_sum_t)
    model_f = model_theta[theta[:,idx]==0]
    model_sum_f = model_sum[theta[:,idx]==0]
    cond_res_f[:,idx] = np.sum(model_f,axis=0)
    cond_res_f_p[:,idx] = (np.sum(model_f,axis=0))/np.sum(model_sum_f)
    
np.save(os.path.join(out_folder, 'cond_res_t.npy'),cond_res_t)
np.save(os.path.join(out_folder, 'cond_res_t_p.npy'),cond_res_t_p)
np.save(os.path.join(out_folder, 'cond_res_f.npy'),cond_res_f)
np.save(os.path.join(out_folder, 'cond_res_f_p.npy'),cond_res_f_p)
