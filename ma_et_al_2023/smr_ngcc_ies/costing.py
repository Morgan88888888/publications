#!/usr/bin/env python
# coding: utf-8

# # today is 0209

# In[1]:


import pandas as pd


# In[2]:


# define the plant configurations
case_list = ['smr_ngsc_ccs','smr_ngcc_ccs','smr_ngsc','smr_ngcc','smr']

class Plant:
    def __init__(self,ptype,pcapacity):
        self.type = ptype
        self.cf = pcapacity
        
    def check(self):
        if self.type not in case_list:
            print("Please select the plant type as one in the list ", case_list)
        else:
            if 0.4 <= self.cf <= 1:
                 print("The plant configurations have been correctly settled")
            else:
                 print("Please select the plant capacity in the range of 40 % to 100 %")
    def files(self):
        title = case_list.index(self.type)
        if title == 0:
            lookf = {
           "smr_ngsc_ccs" :'case1t',
           "ccs_ccs" : 'case1c'   
            }
        elif title == 1:
            lookf = {
            "smr_ngcc_ccs" :'case2t',
            "ccs_ccs" : 'case2c'   
            }
        else:
            lookf = {
            "smr_ngsc" :'case3t',
            "smr_ngcc" : 'case4t',
            'smr': 'smr'
            }
        return lookf
                          


# In[3]:


def read_data(fname):
    df = pd.read_csv(f"result/{fname}.csv", header= 0, index_col = 0)
    data = []
    for cl in df.columns:
        data.append(df[cl].tolist())
    return data


# In[4]:


#Different scaling rules
def rl0(b, cf):
    t = b * cf ** 0.8
    return t

def rl1(b ,cf1, cf2):
    t = b * 0.4 *cf1 ** 0.6 + b *0.6 * cf2 ** 0.6
    return t

def rl2(b, cf):
    t = b * cf ** 0.41
    return t

def rl3(b, cf):
    t = b * cf ** 0.83
    return t

def rl4(b, cf):
    t = b * cf ** 0.79
    return t

def rlv(b, cf):
    t = b * cf 
    return t
    


# In[46]:


def LCOH_cal(summary, product,cf, case = 'smr',co2= 0):
    LCOH = pd.DataFrame()
    LCOH['Value, $/kg H2'] = (summary.iloc[:,0]/(product.iloc[0,1]*365*24/1000000*cf)).round(2)
    if case != 'smr':
        #print(product)
        LCOH.at['CO2 S&T','Value, $/kg H2'] = (co2/1000*24*365/100000/(product.iloc[0,1]*365*24/1000000)).round(2)
    else:
        LCOH.at['CO2 S&T','Value, $/kg H2'] = (194054/0.97*0.9/1000*24*365/100000/(product.iloc[0,1]*365*24/1000000*0.9)).round(2)
    LCOH.at['Total','Value, $/kg H2'] = LCOH['Value, $/kg H2'].sum()
    LCOH.to_csv('LCOH.csv')
    return LCOH


# In[72]:


def cost_cal(product, case, ratio = None):
    cost = pd.read_csv("reference/cost.csv", header= 0, index_col = 0).fillna(0)
    idxl = cost.index

    if case is not None:
        cost[case] = pd.Series(0)
        cost.at[idxl[0], case] = sum(rl0(cost.at[idxl[0], cl], 
                                     product.replace(0,1).at[product.index[case.split('_').index(cl)], case]/product.replace(0,1).at[product.index[case.split('_').index(cl)], cl]) 
                                  for cl in case.split('_') )
        cost.at[idxl[1], case] = sum(rl1(cost.at[idxl[1], cl], 
                                     product.at['fg',case]/product.replace(0,1).at['fg',cl],
                                     product.at['co2',case]/product.replace(0,1).at['co2',cl])                              
                                  for cl in case.split('_'))
        cost.at[idxl[2], case] = sum(rl2(cost.at[idxl[2], cl], 
                                     product.at['co2',case]/product.replace(0,1).at['co2',cl])                              
                                  for cl in case.split('_'))
        cost.at[idxl[3], case] = sum(rl3(cost.at[idxl[3], cl], 
                                     product.at['co2',case]/product.replace(0,1).at['co2',cl])                              
                                  for cl in case.split('_'))
        cost.at[idxl[4], case] = sum(rl4(cost.at[idxl[4], cl], 
                                     product.at['co2',case]/product.replace(0,1).at['co2',cl])                              
                                  for cl in case.split('_'))
        if case == 'smr_ngsc':
            cost.at[idxl[5], case] = (3.87*1.25 +
                                     ((0.0076+0.02/1.25)*(cost.at[idxl[1], case]+cost.at[idxl[2], case]+cost.at[idxl[3], case]+cost.at[idxl[4], case]+480.3131)
                                     )*1.25)
        elif case == 'smr_ngcc':
            cost.at[idxl[5], case] = (4.53 *1.25+
                                     ((0.0076+0.02/1.25)*(cost.at[idxl[1], case]+cost.at[idxl[2], case]+cost.at[idxl[3], case]+cost.at[idxl[4], case]+677.45)
                                     )*1.25)
        else:
            cost.at[idxl[5], case] = sum(rl0(cost.at[idxl[5], cl], 
                                     product.at[product.index[case.split('_').index(cl)], case]/product.replace(0,1).at[product.index[case.split('_').index(cl)], cl]) 
                                  for cl in case.split('_') ) 
        cost.at[idxl[6], case] = sum(rlv(cost.at[idxl[6], cl], 
                                     product.at[product.index[case.split('_').index(cl)], case]/product.replace(0,1).at[product.index[case.split('_').index(cl)], cl]) 
                                  for cl in case.split('_') ) 
        cost.at[idxl[7], case] = sum(rlv(cost.at[idxl[7], cl], 
                                     product.at[product.index[case.split('_').index(cl)], case]/product.replace(0,1).at[product.index[case.split('_').index(cl)], cl]) 
                                  for cl in case.split('_') ) 
        if case in ["smr_ngsc",'smr_ngcc']:
            cost.at[idxl[8], case] = 0
        else:
            cost.at[idxl[8], case] = sum(rlv(cost.at[idxl[8], cl], 
                                     product.at[product.index[case.split('_').index(cl)], case]/product.replace(0,1).at[product.index[case.split('_').index(cl)], cl]) 
                                  for cl in case.split('_') ) 
        if ratio != None:
            cost.at[idxl[8], case] = cost.at[idxl[8], case] * ratio
            
        cost.at[idxl[9], case] = sum(rlv(cost.at[idxl[9], cl], 
                                     product.at['co2',case]/product.replace(0,1).at['co2',cl])                              
                                  for cl in case.split('_'))
        cost.at[idxl[10], case] = sum(rlv(cost.at[idxl[10], cl], 
                                     product.at[product.index[case.split('_').index(cl)], case]/product.replace(0,1).at[product.index[case.split('_').index(cl)], cl]) 
                                  for cl in case.split('_') ) 
    for col in cost.columns:
        cost.at['TOC', col] = sum(cost.at[idx,col] for idx in idxl if "OC" in idx)
    cost = pd.concat([cost,
                      pd.DataFrame((cost.loc["TOC",:] * 1.069802).rename('TASC', inplace = True)).T],
                      axis = 0)
    if case == 'smr':
        cost = pd.concat([cost,
                          pd.DataFrame((cost.loc["TASC",:] * 0.0586).rename('ACC', inplace = True)).T],
                          axis = 0)  
    else:
        cost = pd.concat([cost,
                          pd.DataFrame((cost.loc["TASC",:] * 0.0707).rename('ACC', inplace = True)).T],
                          axis = 0)
    
    cost = pd.concat([cost,
                      pd.DataFrame((cost.loc["ACC",:]+ cost.loc["FOM",:]).rename('TFC', inplace = True)).T],
                      axis = 0)
    for col in cost.columns:
        cost.at['TVC', col] = sum(cost.at[idx,col] for idx in idxl if "VC" in idx)

    return cost


# In[73]:


def cost_sum(product,total,cost,case,cf, i = 0):
    summary = pd.DataFrame()
    if 'ccs' in case:
        case0 = case
        case = 'ccs_ccs'
    else:
        case0 = case
    if 'ccs' not in case0:
        summary.at['Capital',case+str(i)] = cost.at['ACC',case]
        summary.at['Fixed O&M',case+str(i)] = - cost.at['ACC',case] + cost.at['TFC',case]
        summary.at['Variable O&M',case+str(i)] = cost.at['TVC',case]*cf
        summary.at['Fuel',case+str(i)] = -total.at['ng',case]*cf
    else:
        summary.at['FC',case0+str(i)] = -cost.at['TFC',case]/2 \
                                 -sum(rl0(cost.at['TFC', cl], 
                                     product.replace(0,1).at[product.index[case0.split('_').index(cl)], case0]/product.replace(0,1).at[product.index[case0.split('_').index(cl)], cl]) 
                                      for cl in case0.split('_')[:-1])
        summary.at['VC',case0+str(i)] = -sum(rlv(cost.at['TVC', cl], 
                                     product.at[product.index[case0.split('_').index(cl)], case0]/product.replace(0,1).at[product.index[case0.split('_').index(cl)], cl]) 
                                     for cl in case0.split('_')[:-1])  \
                                -cost.at['TVC',case]/2 \
                                + total.at['ng',case0]
    return summary


# In[74]:


def costing(plant):

    co2list = {'smr_ngsc' : (123947.67+98925.46),
		   'smr_ngcc' : (123956.43+181218.76)}
    
    case = plant.type
    cf = plant.cf
    if 'ccs' in case:
        case0 = case
        case = 'ccs_ccs'
    else:
        case0 = case

    for i in range(len(read_data(plant.files()[case0]))):    
        product = pd.read_csv("reference/product.csv", header= 0, index_col = [0,1])
        if case is not None: 
            product[case] = read_data(plant.files()[case])[i]
        product.loc[['co2']] = product.loc[['co2']] /1000
        product.loc[['ng']] = product.loc[['ng']] * 49565.84/1e6
        product.reset_index(inplace = True)
        product.set_index('item',inplace = True,drop = True)
        price = pd.read_csv('reference/price.csv', header = 0, index_col = 0)
        total = pd.DataFrame()
        for col in product.drop('unit', axis = 1).columns:
            for index in price.instant.index:
                total.at[index,col] = product[col][index] * price.instant[index]*24*365/1e6
	   
        cost= cost_cal(product,case)

        if 'ccs' in case0:
            product = pd.read_csv("reference/product.csv", header= 0, index_col = [0,1])
            product[case0] = read_data(plant.files()[case0])[i]
            product.loc[['co2']] = product.loc[['co2']] /1000
            product.loc[['ng']] = product.loc[['ng']] * 49565.84/1e6
            product.reset_index(inplace = True)
            product.set_index('item',inplace = True,drop = True)
            total = pd.DataFrame()
            for col in product.drop('unit', axis = 1).columns:
                for index in price.instant.index:
                    total.at[index,col] = product[col][index] * price.instant[index]*24*365/1e6
	    #var and fixed cost

    summary = cost_sum(product,total,cost,case,cf)
    try:
        co2 = co2list[case]
    except:
        co2 = 0
    LCOH = LCOH_cal(summary, product,cf,case, co2)   
    #print(LCOH)
    return LCOH

plant = Plant('smr_ngsc', 1)




# In[ ]:




