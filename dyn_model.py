# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:14:48 2019

@author: cwmorin & eburnor
"""

import numpy as np
import matplotlib.pyplot as plt
from detectProb import detectProbMC
import time
import math as math
import random as random
import pathlib as pathlib

divider = 10

class branch:
    # initialize static branch parameters
    def __init__(self,idnum, vel, leng, pop, prev, parents, resolution):
        
        self.id = int(idnum)
        self.res = resolution
     
        
        self.length_meters = leng # meters
        
        segments_res = leng/self.res
        if segments_res < 1:
            self.length = 1
        else:
            self.length = round(segments_res)     
                                   
        self.population = round(pop)
        self.prevalence = prev
        
        self.parents = parents
        self.parents = [int(i) for i in self.parents]
        
        self.no_parents = True
        for i in range(0,len(self.parents)):
            if self.parents[i] >= 0:
                self.no_parents = False
        
    # create copies for updating and created histories to retain most recent days data
    def initializeHistory(self,sims):
        self.pathSection = np.zeros((int(self.length)))
        self.pathSectCpy = np.zeros((int(self.length)))
        self.floSection = np.zeros((int(self.length)))
        self.floSectCpy = np.zeros((int(self.length)))
        self.conSection = np.zeros((int(self.length)))

        self.pathHistory = np.zeros((sims,24,int(self.length)))
        self.floHistory = np.zeros((sims,24,int(self.length)))
        self.conHistory = np.zeros((sims,24,int(self.length)))

        
        self.detectMean = np.zeros((24,int(self.length)))
        self.detectSTD = np.zeros((24,int(self.length)))
        self.detectLower = np.zeros((24,int(self.length)))
        self.detectUpper = np.zeros((24,int(self.length)))
        
        
        self.detectMedianPM = np.zeros((24,int(self.length)))
   
        self.detectMeanPM = np.zeros((24,int(self.length)))
        self.detectSTDPM= np.zeros((24,int(self.length)))
        self.detectLowerPM = np.zeros((24,int(self.length)))
        self.detectUpperPM = np.zeros((24,int(self.length)))
        
        self.detectMedian = np.zeros((24,int(self.length)))
        self.detectMedianPM = np.zeros((24,int(self.length)))
        
    # initialize distributed parameters 
    def initializeParams(self, defR, fl, eflow, decR, probShed, shedDist):
        
        # Create array for people along branch 
        self.pop = np.zeros((self.population, 6))
        
        # Randomly assign subset of branch to be infected with Typhoid, based on estimated prevalence
        self.pop[:,0] = random.choices([0,1], weights = (1-self.prevalence, self.prevalence), k = self.population)
    
        # Set up population along branch
        for row in self.pop:
            row[1] = np.random.normal(defR[0],defR[1])

            if row[1] < 0:
                row[1] = 0
            
            if row[0] == 1:
                
                # If the person is sick, randomly choose days since infection with equal weights
                days = random.choices(range(14), k = 1)[0]
                row[2] = days
                probShedPerson = probShed[days, 1]
                row[3] = probShedPerson
                
                shedding = random.choices([0,1], weights = (1-probShedPerson, probShedPerson), k = 1)[0]
                if shedding == 1:
                    row[4]= np.random.choice(shedDist, 1, replace = True)[0]
                row[5] = row[1]*row[4]
                
        #print(self.id)
        self.shedRate = np.sum(self.pop[:,5])/(int(self.length))
        #print(self.shedRate)
        
        self.flow = np.random.normal(fl[0],fl[1]) 
        # flow cannot be negative
        if(self.flow < 0):
            self.flow = 0
     
        self.flowRate = (self.population*self.flow)/(int(self.length)) # Liters / day in the section of the branch
     
        # Liters per hour
        self.eflow = np.random.normal(eflow[0], eflow[1])
        if(self.eflow < 0):
            self.eflow = 0
            
        self.velocity = ((self.eflow + self.flowRate/24) * (1/1000) * (1/0.5)) / self.res
      
        self.decRate = decR
    
    #helper function to save most recent day's data
    def copyToFloHist(self,sim,hour):
        self.floHistory[sim,hour,:] = self.floSection
    
    #helper function to save most recent day's data    
    def copyToPathHist(self,sim,hour):
        self.pathHistory[sim,hour,:] = self.pathSection
     
    #helper function to save most recent day's data    
    def copyToConHist(self,sim,hour):
        self.conHistory[sim,hour,:] = self.conSection
    
    #calculate statistics across simulations and use them to calculate detection statistics
    def calcStatistics(self,samp,sens,bact, sims, conc_only):
        
        if(conc_only):
            
            self.simMean = np.mean(self.conHistory,axis=0)
            self.simMedian = np.median(self.conHistory,axis=0)
            self.simSTD = np.std(self.conHistory,axis=0)
            self.simSE = self.simSTD / math.sqrt(sims)
            self.simMedian = np.median(self.conHistory,axis=0)
            
        else:
            self.simMean = np.mean(self.conHistory,axis=0)    
            self.simSTD = np.std(self.conHistory,axis=0)
            self.simSE = self.simSTD / math.sqrt(sims)
            self.simMedian = np.median(self.conHistory,axis=0)
    
            for i in range(24):
               for j in range(int(self.length)):
               
                   # Temp: 
                   detection = detectProbMC(self.simMean[i,j],samp,sens,bact)
              
                   self.detectMean[i,j] = detection[0]
                   self.detectSTD[i,j] = detection[1]
                   self.detectLower[i,j] = detection[2]
                   self.detectUpper[i,j] = detection[3]
                   self.detectMedian[i,j] = detection[4]
                   
                   
     # Calc stats using concentration median
    def calcStatisticsMedian(self,samp,sens,bact, sims, conc_only):
        
        if(conc_only):
             self.simMedian = np.median(self.conHistory,axis=0)
             self.simMean = np.mean(self.conHistory,axis=0)    
             self.simSTD = np.std(self.conHistory,axis=0)
             self.simSE = self.simSTD / math.sqrt(sims)
             
        else:
             
             self.simMedian = np.median(self.conHistory,axis=0)
             self.simMean = np.mean(self.conHistory,axis=0)    
             self.simSTD = np.std(self.conHistory,axis=0)
             self.simSE = self.simSTD / math.sqrt(sims)
     
             for i in range(24):
                for j in range(int(self.length)):
    
                    detection = detectProbMC(self.simMedian[i,j],samp,sens,bact)
               
                    self.detectMean[i,j] = detection[0]
                    self.detectSTD[i,j] = detection[1]
                    self.detectLower[i,j] = detection[2]
                    self.detectUpper[i,j] = detection[3]
                    self.detectMedian[i,j] = detection[4]
                   
          
             

#ws = dyn.WaterSystem(fixParams=fixed_p, distParams = dist_p, dayParams = day_p, probShed = probShed, sheddingDist = shedDist, samVol = 1, sens = 0.83, bactNum = 10)
class WaterSystem:
    
    def __init__(self,fixParams,distParams,dayParams, probShed, sheddingDist, samVol, sens, bactNum, savePath, resolution):
        
        self.filePath = savePath
        self.fixed = fixParams
        self.resolution = resolution
        num_rows, num_cols = self.fixed.shape
        
        dist_p = distParams
        day_p = dayParams
        
        self.flowCycle = day_p[:,0]
        self.pathCycle = day_p[:,1]
        
        # Set up detection parameters
        self.sampleVol = samVol
        self.sensitivity = sens
        self.bacteriaNum = bactNum 
        
        # Set up system distributed paramters
        self.def_r = (dist_p[0],dist_p[1])
        self.fl = (dist_p[2],dist_p[3])
        self.efl = (dist_p[7], dist_p[8])
        self.des_r = (dist_p[4],dist_p[5],dist_p[6])
        self.probShed = probShed
        self.shedDist = sheddingDist

        # determine number of branches in file
        self.numBranches = len(self.fixed[:,0])
        
        # detemine length of longest branch (in sections)
        self.longestBranch = int(np.max(self.fixed[:,2])/resolution)
        
        # Branches needs to be a dictionary
        #create dictionary 
        self.branches= {}
        
        #loop to add all branches to list
        for i in range(self.numBranches):
            
            branch_id = int(self.fixed[i][0])
            
            parents = self.fixed[i,5:num_cols]
            parents = parents.tolist()

            self.branches[branch_id] = branch(branch_id, self.fixed[i][1], self.fixed[i][2], 
                                        self.fixed[i][3], self.fixed[i][4], parents, resolution)

      
       
    def pathogenLoad(self,time, i, bran, hour):
   
       dist = int(round(time*bran.velocity))
       
       #print(bran.id, " ", time, " ", i, " ", hour, " ", dist)
       
       if(dist <= i):
           #print(str(bran.id) + " " + str(i) + " " + str(dist))
           path = bran.pathSection[i-dist]*bran.decRate + time*self.pathCycle[hour]*bran.shedRate
           return path
       elif (bran.no_parents == True):
           newTime = i / (bran.velocity)
           path = newTime*self.pathCycle[hour]*bran.shedRate
           return path
       else:
           timeBranch = i / (bran.velocity)
           timeExtra = time - timeBranch
           path0 = timeBranch*self.pathCycle[hour]*bran.shedRate
               
           parent_path_total = 0
           for k in range(0, len(bran.parents)):
               parent_id = bran.parents[k]
               if(parent_id == -1):
                   parent_path_total += 0
               else:
                   parent_branch = self.branches[parent_id]
                   path = self.pathogenLoad(timeExtra,int(parent_branch.length)-1,parent_branch,hour)
                   parent_path_total += path
           path = path0 + parent_path_total
           return path
        
     
      
    def floLoad(self,time, i, bran, hour):
        
        dist = int(round(time*bran.velocity)) # distance (in terms of sections)
        #print(bran.id, " ", time, " ", i, " ", hour, " ", dist)
        if(dist <= i):
            flo = bran.floSection[i-dist] + time*self.flowCycle[hour]*bran.flowRate + bran.eflow
            return flo
        elif (bran.no_parents == True):
            newTime = i / (bran.velocity)
            flo = newTime*self.flowCycle[hour]*bran.flowRate + newTime*bran.eflow
            return flo
        else:
            timeBranch = i / (bran.velocity)
            timeExtra = time - timeBranch
            
            flo0 = timeBranch*self.flowCycle[hour]*bran.flowRate + timeBranch*bran.eflow
            
            parent_flo_total = 0
            for k in range(0, len(bran.parents)):
                parent_id = bran.parents[k]
                if(parent_id == -1):
                    parent_flo_total += 0
                else:
                    parent_branch = self.branches[parent_id]
                    flo = self.floLoad(timeExtra,int(parent_branch.length)-1,parent_branch,hour)
                    parent_flo_total += flo
            flo_total = flo0 + parent_flo_total
            return flo_total
       

    def update(self, hour):
        for bran in self.branches.values():
            
            for i in range(len(bran.pathSection)-1,-1, -1):
                
                bran.pathSectCpy[i] = self.pathogenLoad(1,i,bran,hour)
                bran.floSectCpy[i] = self.floLoad(1,i,bran,hour)
                
                #print(self.pathogenLoad(1,i,bran,hour))
                #print(self.floLoad(1,i,bran,hour))
                
                if bran.floSectCpy[i] == 0:
                    bran.conSection[i] = 0
                else:
                    bran.conSection[i] = bran.pathSectCpy[i]/bran.floSectCpy[i]
        for bran in self.branches.values():
            np.copyto(bran.pathSection,bran.pathSectCpy)
            np.copyto(bran.floSection,bran.floSectCpy)
            

                    
    def run(self,days,sims):
        
        self.simulations = sims
        # Initialize all branches
        for bran in self.branches.values():
            bran.initializeHistory(sims)
        # Loop through each simulation

        for i in range(sims):
            start1 = time.time()
       
            self.desiccation_r = random.choices([self.des_r[0],self.des_r[1],self.des_r[2]], k = 1)[0]
            

          
            for bran in self.branches.values():
                    bran.initializeParams(self.def_r, self.fl, self.efl, self.desiccation_r, self.probShed, self.shedDist)
        
            # Run each simulation over specified number of days
            
            for j in range(days):
                for k in range(24):
                    # Update model hourly
                    self.update(k)
                    # Copy updates to history 
                    for bran in self.branches.values():
                        bran.copyToPathHist(i,k)
                        bran.copyToFloHist(i,k)
                        bran.copyToConHist(i,k)
            
            stop1 = time.time()
            print("one simulation", stop1-start1)
        
        calc_start = time.time()
        # When all simulations are complete, compute statistics 
        for bran in self.branches.values():
            bran.calcStatisticsMedian(self.sampleVol,self.sensitivity,self.bacteriaNum, sims, False)
        calc_end = time.time()
        print("Calc Statistics ", calc_end-calc_start)
    
    def reComputeStatistics(self, sampleVol, sensitivity, bacteriaNum, sims): # Recompute statistics using concentration mean
    
        for bran in self.branches.values():
            bran.calcStatistics(sampleVol,sensitivity,bacteriaNum, sims, False)
            
    def reComputeStatisticsMedian(self, sampleVol, sensitivity, bacteriaNum, sims): # Recompute statistics using concentration median
     
         for bran in self.branches.values():
             bran.calcStatisticsMedian(sampleVol,sensitivity,bacteriaNum, sims, False)
            
    def dayTest(self):
        
        for bran in self.branches.values():
            bran.initializeHistory(1)
        # Loop through each simulatio
        
        self.desiccation_r = random.choices([self.des_r[0],self.des_r[1],self.des_r[2]], k = 1)[0]
        
        for bran in self.branches.values():
            bran.initializeParams(self.def_r, self.fl, self.efl, self.desiccation_r, self.probShed, self.shedDist)
            
        
        history_rows = 0
        for bran in self.branches.values():
            sections = bran.length
            print(sections)
            history_rows += sections
            
        prevDayConc = np.zeros((history_rows, 24))
        
        print("array dimensions ", np.shape(prevDayConc))
        
        print("Start Day 1")
        day = 0
        mad = 1
        while(mad > 0):
            day += 1
            for k in range(24):
                # Update model hourly
                self.update(k)
                # Copy updates to history 
                for bran in self.branches.values():
                    bran.copyToPathHist(0,k)
                    bran.copyToFloHist(0,k)
                    bran.copyToConHist(0,k)
             
        
    
            # When all simulations are complete, compute statistics 
            for bran in self.branches.values():
                bran.calcStatistics(self.sampleVol,self.sensitivity,self.bacteriaNum, 1, True)
            
            i = 0
            for bran in self.branches.values():
                idnum = np.array([bran.id]*(int(bran.length)))
                seg = np.array(range(0,int(bran.length)))
                

                
                if(i == 0):
                    concatCmean = np.column_stack((idnum,seg,np.transpose(bran.simMean)))
                else:
                    branCMean = np.column_stack((idnum,seg,np.transpose(bran.simMean)))
                    concatCmean = np.concatenate((concatCmean, branCMean), axis = 0)
                i += 1
            
            todayConc = concatCmean[:,2:26]
            diff = np.subtract(todayConc, prevDayConc)
            absDiff = np.absolute(diff)
            
            mad = np.sum(absDiff)/np.size(absDiff)
            prevDayConc = todayConc
            print("Absolute Difference ", mad)
        print(day) 
        return(concatCmean)
            
            
    def graphByHour(self, saveFile):
         minMax = self.findMinMaxF()
         fig = plt.figure(figsize=(7,8))
         for i in range(24):
             values = np.empty((int(self.longestBranch),int(self.numBranches)))
             values[:] = np.NaN
             j = 0
             for bran in self.branches.values():
                 values[0:bran.length,j] = bran.detectMean[i,:]
                 j = j + 1
             ax = fig.add_subplot(6,4,i+1)
             im = ax.pcolormesh(values, cmap='hot_r',vmin=minMax[0],vmax=minMax[1])
             ax.set_facecolor((.7,.7,.7))
             if((i+1)<21):
                 plt.xticks([])
             if(((i+1)%4) != 1):
                 plt.yticks([])
             r = int(i + 1)
             ax.set_title('Hour = %r' %r)
         cbaxes = fig.add_axes([.15, 0.07, 0.7, 0.025])
         fig.colorbar(im,cax=cbaxes,orientation="horizontal")
         fig.subplots_adjust(left=.04,bottom=.13,right=.98,top=.95,wspace=0.04,hspace=0.3)
         
         if(saveFile):
             savePath = pathlib.Path(self.filePath, "PlotByHour.tiff")
             plt.savefig(savePath)
    
         else:
            return(fig)
     
    def graphByBranch(self, saveFile):
        
         
         minMax = self.findMinMaxF()
         
         if(saveFile):
             for bran in self.branches.values():
                 hours = np.arange(24)
                 sections = np.arange(int(bran.length))
                 x,y = np.meshgrid(hours,sections)
                 #minVal = np.min(bran.detectMean)
                 #maxVal = np.max(bran.detectMean)
                 values = np.transpose(bran.detectMean)
                 fig,ax = plt.subplots()
                 im = ax.pcolormesh(values, cmap='hot_r',vmin=minMax[0],vmax=minMax[1])
                 fig.colorbar(im, ax=ax)
                 plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
                 brID = 'Branch ' + str(bran.id)
                 ax.set_title(brID)
                 ax.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
                 ax.set_xlabel('Hour')
                 ax.set_ylabel('Section (' + str(self.resolution) + ' m)')
                 plt.tight_layout()
                 filename = brID + 'detect.tiff'
                 savePath = pathlib.Path(self.filePath, filename)
                 plt.savefig(savePath)
                 plt.close(fig)
         else:
            key = list(self.branches.keys())[1]
            bran = self.branches[key]
            hours = np.arange(24)
            sections = np.arange(int(bran.length))
            x,y = np.meshgrid(hours,sections)
            #minVal = np.min(bran.detectMean)
            #maxVal = np.max(bran.detectMean)
            values = np.transpose(bran.detectMean)
            fig,ax = plt.subplots()
            im = ax.pcolormesh(values, cmap='hot_r',vmin=minMax[0],vmax=minMax[1])
            fig.colorbar(im, ax=ax)
            plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
            brID = 'Branch ' + str(bran.id)
            ax.set_title(brID)
            ax.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
            ax.set_xlabel('Hour')
            ax.set_ylabel('Section (' + str(self.resolution) + ' m)')
            plt.tight_layout()
            return(fig)
             
         
             
             
    def graphByBranchConc(self, saveFile):
         minMax = self.findMinMaxConc()
         
         if(saveFile):
             for bran in self.branches.values():
                 hours = np.arange(24)
                 sections = np.arange(int(bran.length))
                 x,y = np.meshgrid(hours,sections)
                 #minVal = np.min(bran.simMean)
                 #maxVal = np.max(bran.simMean)
                 values = np.transpose(bran.simMean)
                 #values = np.log10(values)
                 fig,ax = plt.subplots()
                 im = ax.pcolormesh(values, cmap='hot_r',vmin=minMax[0],vmax=minMax[1])
                 fig.colorbar(im, ax=ax)
                 plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
                 brID = 'Branch_' + str(bran.id)
                 ax.set_title(brID)
                 ax.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
                 ax.set_xlabel('Hour')
                 ax.set_ylabel('Section (' + str(self.resolution) + ' m)')
                 plt.tight_layout()
            
                 filename = brID + 'conc.tiff'
                 savePath = pathlib.Path(self.filePath, filename)
                 plt.savefig(savePath)
                 plt.close(fig)
         else:
             key = list(self.branches.keys())[1]
             bran = self.branches[key]
             hours = np.arange(24)
             sections = np.arange(int(bran.length))
             x,y = np.meshgrid(hours,sections)
             #minVal = np.min(bran.simMean)
             #maxVal = np.max(bran.simMean)
             values = np.transpose(bran.simMean)
             values = np.log10(values)
             fig,ax = plt.subplots()
             im = ax.pcolormesh(values, cmap='hot_r',vmin=minMax[0],vmax=minMax[1])
             fig.colorbar(im, ax=ax)
             plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
             brID = 'Branch_' + str(bran.id)
             ax.set_title(brID)
             ax.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
             ax.set_xlabel('Hour')
             ax.set_ylabel('Section (' + str(self.resolution) +  ' m)')
             plt.tight_layout()
             return(fig)
        
            
     
     
    def graphDetectByBranchSection(self,branchID,section,saveFile):
         bran = self.branches[branchID]
         sectM = bran.detectMean
         # sectSD = bran.detectSTD
         sectLower = bran.detectLower
         sectHigher = bran.detectUpper
         
         sectMH = sectM[:,section]
         #sectSDH = sectSD[:,section]
         sectLH = sectLower[:,section]
         sectHH = sectHigher[:,section]

         hour = np.arange(24)
         #plt.errorbar(hour,sectMH, yerr = (ybot, ytop))
         fig = plt.plot(hour, sectMH, '-')
         plt.fill_between(hour, sectLH, sectHH, alpha=0.2)
         plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
         brID = 'Branch ' + str(branchID)+' Section '+str(section)
         plt.title(brID)
         plt.xlabel('Hour')
         plt.ylabel('P(detection)')
         plt.tight_layout()
         
         if(saveFile):
             filename = brID + '.tiff'
             savePath = pathlib.Path(self.filePath, filename)
             plt.savefig(savePath)
             plt.close()
         else:
             return(fig)
         
    def graphConcByBranchSection(self,branchID,section, saveFile):
         bran = self.branches[branchID]
         sectM = bran.simMean
         # sectSD = bran.simSTD
         sectSE = bran.simSTD / math.sqrt(self.simulations)
         
         sectMH = sectM[:,section]
         sectSEH = sectSE[:,section]
         
         ytop= sectMH + 1.96*sectSEH
         ybot = sectMH - 1.96*sectSEH
         for i in range(len(ybot)):
             if ybot[i] < 0:
                 #ybot[i] = abs(0-sectMH[i])
                 ybot[i] = 0
                 
         hour = np.arange(24)
         
         fig = plt.plot(hour, sectMH, '-')
         plt.fill_between(hour, ybot, ytop, alpha=0.2)
         # plt.errorbar(hour,sectMH,yerr = (ybot, ytop))
         plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
         brID = 'Branch ' + str(bran.id)+' Section '+str(section)
         plt.title(brID)
         plt.xlabel('Hour')
         plt.ylabel('Concentration (CFU/L)')
         plt.tight_layout()
         
         if(saveFile):
             filename = brID + '.tiff'
             savePath = pathlib.Path(self.filePath, filename)
             plt.savefig(savePath)
         else:
             return(fig)
         
     

    def findMinMaxF(self):
         minV = 0
         maxV = 0
         for bran in self.branches.values():
             if(np.max(bran.detectMean)>maxV):
                 maxV = np.max(bran.detectMean)
             if(np.min(bran.detectMean)<minV):
                 minV = np.min(bran.detectMean)
         return (minV,maxV)
     
    def findMinMaxConc(self):
         minC = 0
         maxC = 0
         for bran in self.branches.values():
             if(np.max(bran.simMean)>maxC):
                 maxC = np.max(bran.simMean)
             if(np.min(bran.simMean)<minC):
                 minC = np.min(bran.simMean)
         return (minC,maxC)
     
    def createCSV(self,fName):
         #hdr = 'Branch,Section,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,H15,H16,H17,H18,H19,H20,H21,H22,H23,H24'
         
         #writePath = pathlib.Path(self.filePath, filename)
         with open((self.filePath+'Cmean.csv'),'w') as f1,\
         open((self.filePath+'Cstd.csv'),'w') as f2,\
         open((self.filePath+'Dmean.csv'),'w') as f3,\
         open((self.filePath+'Dlower.csv'),'w') as f4,\
         open((self.filePath+'Dupper.csv'),'w') as f5,\
         open((self.filePath+'Cmedian.csv'),'w') as f6,\
         open((self.filePath+'Dmedian.csv'),'w') as f7:
             i = 0
             for bran in self.branches.values():
                 idnum = np.array([bran.id]*(int(bran.length)))
                 seg = np.array(range(0,int(bran.length)))
                 concatCmean = np.column_stack((idnum,seg,np.transpose(bran.simMean)))
                 concatCstd = np.column_stack((idnum,seg,np.transpose(bran.simSTD)))
                 concatDmean = np.column_stack((idnum,seg,np.transpose(bran.detectMean)))
                 concatDlower = np.column_stack((idnum,seg,np.transpose(bran.detectLower)))
                 concatDupper = np.column_stack((idnum,seg,np.transpose(bran.detectUpper)))
                 concatCmedian = np.column_stack((idnum,seg,np.transpose(bran.simMedian)))
                 concatDmedian = np.column_stack((idnum,seg,np.transpose(bran.detectMedian)))
                 
                 if i>0:
                     np.savetxt(f1,concatCmean,delimiter=',',fmt='%.4f')
                     np.savetxt(f2,concatCstd,delimiter=',',fmt='%.4f')
                     np.savetxt(f3,concatDmean,delimiter=',',fmt='%.4f')
                     np.savetxt(f4,concatDlower,delimiter=',',fmt='%.4f')
                     np.savetxt(f5,concatDupper,delimiter=',',fmt='%.4f')
                     np.savetxt(f6,concatCmedian,delimiter=',',fmt='%.4f')
                     np.savetxt(f7,concatDmedian,delimiter=',',fmt='%.4f')
                 else:
                     np.savetxt(f1,concatCmean,delimiter=',',fmt='%.4f')
                     np.savetxt(f2,concatCstd,delimiter=',', fmt='%.4f')
                     np.savetxt(f3,concatDmean,delimiter=',', fmt='%.4f')
                     np.savetxt(f4,concatDlower,delimiter=',', fmt='%.4f')
                     np.savetxt(f5,concatDupper,delimiter=',', fmt='%.4f')
                     np.savetxt(f6,concatCmedian,delimiter=',',fmt='%.4f')
                     np.savetxt(f6,concatCmedian,delimiter=',',fmt='%.4f')
                     np.savetxt(f7,concatDmedian,delimiter=',',fmt='%.4f')
           
# self, fixParams, distParams, dayParams, shedding, samVol, sens, bactNum)

"""

 np.savetxt(f1,concatCmean,delimiter=',',header=hdr,fmt='%.4f')
 np.savetxt(f2,concatCstd,delimiter=',',header=hdr,fmt='%.4f')
 np.savetxt(f3,concatDmean,delimiter=',',header=hdr,fmt='%.4f')
 np.savetxt(f4,concatDlower,delimiter=',',header=hdr,fmt='%.4f')
 np.savetxt(f5,concatDupper,delimiter=',',header=hdr,fmt='%.4f')
path = '/Users/elisabethburnor/Desktop/Typhoid model GIS/sample_runs/simple_run/'

# Read in csvs
fixed_p = np.loadtxt('/Users/elisabethburnor/Desktop/Typhoid model GIS/sample_runs/simple_run/Fixed_Params.csv',skiprows = 1, delimiter=",")
dist_p = np.loadtxt('/Users/elisabethburnor/Desktop/Typhoid model GIS/sample_runs/simple_run/Dist_Params.csv',skiprows = 1, delimiter=",")
day_p = np.loadtxt('/Users/elisabethburnor/Desktop/Typhoid model GIS/sample_runs/simple_run/Day_Params.csv',skiprows = 1, delimiter=",")
sheddingEst = np.loadtxt('/Users/elisabethburnor/Desktop/Typhoid model GIS/sample_runs/simple_run/shedding.csv',skiprows = 1, delimiter=",")

# Set up water system
# WaterSystem(fixParams,distParams,dayParams, shedding, samVol, sens, bactNum)
"""
#ws = WaterSystem(fixParams=fixed_p, distParams = dist_p, dayParams = day_p, shedding = sheddingEst, samVol = 1, sens = 0.83, bactNum = 10)
#ws.dayTest()



#example.graphByBranch()
#example.graphByBranchSection(2,2)

"""
example.graphByHour()
example.graphDetectByBranchSection(5,5)
example.graphConcByBranchSection(5,5)
example.graphDetectByBranchSection(13,4)
example.graphConcByBranchSection(13,4)
example.graphByBranch()
example.graphByBranchConc()
"""




# a.createCSV('run1')
# a.graphByHour()
#example.graphDetectByBranchSection(5,5)
#example.graphConcByBranchSection(5,5)
#example.graphDetectByBranchSection(13,4)
#example.graphConcByBranchSection(13,4)
#example.graphByBranch()
#example.graphByBranchConc()

#for bran in a.branches.values():
#    print(bran.length)

# fixed = a.fixed
