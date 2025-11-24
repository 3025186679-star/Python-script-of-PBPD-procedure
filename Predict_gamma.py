import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
import numpy as np
import joblib
import pickle
from scipy.stats import lognorm
#Define PDF function
def pdf(x_gm,ls,lm):
    return (1/(x_gm*ls*np.sqrt(2*3.14159)))*(np.exp(-((np.log(x_gm)-lm)**2)/(2*(ls)**2)))  
#Calculate the inverse x for a given probability level p
def find_x_for_probability(p, ls, lm):
    return lognorm.ppf(p, s=ls, scale=np.exp(lm))
def predict(a1,a2,a3,ξ,ξ1,ξ2,μ,T,p):
    ###########################################
    ##Energy modification coefficient under MCE
    ###########################################
    #Load XGBoost model
    with open("XGBoost_model\\best_xgb_model_biaozhuncha.pkl", 'rb') as f:
        model_biaozhuncha = pickle.load(f)
    with open("XGBoost_model\\best_xgb_model_junzhi.pkl", 'rb') as f:
        model_junzhi = pickle.load(f)
    scaler = joblib.load("XGBoost_model\\scaler_biaozhuncha.pkl")###Load the saved scaler
    inputs =[a1,a2,a3,ξ,ξ1,ξ2,μ,T]
    #Convert input to DataFrame
    input_data = pd.DataFrame([inputs], columns=['α1', 'α2', 'α3', 'ξ', 'ξ1', 'ξ2', 'μ','T'])
    #Standardize input data
    input_data_scaled = scaler.transform(input_data)
    #Predict mean and standard deviation
    prediction_biaozhuncha = model_biaozhuncha.predict(input_data_scaled)
    prediction_junzhi = model_junzhi.predict(input_data_scaled)
    #Predict Energy modification coefficient under MCE
    gama_MCE_pre=find_x_for_probability(p,prediction_biaozhuncha,prediction_junzhi)
    ###########################################
    ##Energy modification coefficient under DBE
    ##Reference: X. Zhou, Y. Chen, K. Ke, M.C.H. Yam, H. Li, Hybrid steel staggered truss frame (SSTF): A probabilistic spectral energy modification coefficient surface model for damage-control evaluation and performance insights, Journal of Building Engineering 45 (2022) 103556. https://doi.org/10.1016/j.jobe.2021.103556.
    ###########################################
    #Calculate fitting formula parameters
    a_DBE = np.exp(-7.424/ξ1 - 3.051*a1 + 1.065)
    b_DBE = 1/(-0.609*np.sqrt(ξ1) + 0.758*np.sqrt(a1) - 0.363)
    c_DBE = 1/(-1.506/np.log(ξ1) + 0.859/np.log(a1) + 0.526)
    d_DBE = 0.028/ξ1 + 0.025*a1**1.5 - 0.024
    e_DBE = -0.067/ξ1 - 0.090*a1**1.5 + 0.072
    f_DBE = np.exp(-2.094/(ξ1**1.5) - 1.258*a1**1.5 - 0.168)
    #Calculate mean and standard deviation
    log_mean_DBE = a_DBE*T**b_DBE + c_DBE
    log_std_DBE = d_DBE*T**2 + e_DBE*T + f_DBE
    #Generate data points
    xxxx = np.arange(0.0001, 10.001, 0.001)  #Energy modification coefficient range
    #Calculate Probability Density Function (PDF)
    PDF_x_DBE = (1/xxxx) * (1/(np.sqrt(2*np.pi)*log_std_DBE)) * \
            np.exp(-((np.log(xxxx) - log_mean_DBE)**2 / (2*log_std_DBE**2)))
    #Calculate Cumulative Distribution Function (CDF)
    CDF_x_DBE = cumulative_trapezoid(PDF_x_DBE, xxxx, initial=0)
    #Interpolate to find the energy coefficient corresponding to the specified probability
    gama_DBE_pre = np.interp(p, CDF_x_DBE, xxxx)
    #Return the Energy modification coefficient
    return(gama_DBE_pre,gama_MCE_pre)