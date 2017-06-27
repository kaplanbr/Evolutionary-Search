# -*- coding: utf-8 -*-
"""
@author: KBK
"""
import numpy as np
import time
from sklearn.model_selection import cross_val_score
import pandas as pd
                                     
class EvolutionaryFeatureSubsetting(object):

    def __init__(self,genelength,populationsize,composition=[.2,.2,.1,.1,.4]):
        self.genelength = genelength 
        self.populationsize = populationsize 
        self.population = {}
        self.nbestfits = int(populationsize*composition[0])
        self.nbestcross = int(populationsize*composition[1])
        self.nbestmutants = int(populationsize*composition[2])
        self.nbestcrossmutant = int(populationsize*composition[3])
        self.nnewborns = int(populationsize*composition[4])
        self.manid = 0
        self.generation = 0
        self.bestfitofgen = {}
        
    def initiateLife(self,n=None):
        g = self.generation
        if self.generation==0:
            n = self.populationsize
            self.population[self.generation] = {}
        for man in range(g*n,g*n+n):
            gene = [np.random.binomial(1,0.5) for i in range(self.genelength)]
            self.population[self.generation][self.manid] = [gene,0]
            self.manid+=1
        print "___________Generation", g ,"____________"
        print n, "random livings are born"
        
    def nextGeneration(self,population):
        tmp_dict = self.population[self.generation]
        self.ranking = sorted([[m,tmp_dict[m]] for m in tmp_dict.keys()],key=lambda x:x[1][1],reverse=True)[0:self.nbestfits]
        self.bestfitofgen[self.generation] = [self.ranking[0][0],self.ranking[0][1][1]]
        print "Best score was", self.ranking[0][1][1],"of id:",self.ranking[0][0]
        self.bestfits = [m[0] for m in self.ranking]
        self.generation = population.keys()[-1]+1
        self.population[self.generation] = {}
        self.initiateLife(n=self.nnewborns)
        for bf in self.bestfits:
            self.population[self.generation][bf] = self.population[self.generation-1][bf]    
        print "%i survived from previous geneneration" %self.nbestfits           
        self.bestcrosses = self.crossover(self.bestfits,returnids=True)
        print "%i bestfit crossovers are born" %self.nbestcross           
        self.mutate(self.bestfits,0.05,bestfitmutants=True)
        print "%i mutated bestfits are born" %self.nbestmutants           
        self.mutate(self.bestcrosses,0.05)
        print "%i mutated bestfit crossovers are born" %self.nbestcrossmutant               
    
     
    def fitnessScore(self,X,y,estimator,search_estimator=False):
        if search_estimator:
            estimator.fit(X,y)
            score = estimator.best_score_
        else:
            score = np.mean(cross_val_score(estimator,X, y, cv = 3, scoring = "r2"))
        return score
      
    def fit(self,X,y,estimator,ngenerations,search_estimator=False):
        print "%i generations of evolution has started" %ngenerations
        self.initiateLife()
        n = ngenerations
        for gen in range(n):
            t = time.time()
            g = self.generation
            for man in self.population[g].keys():
                chosen_columns = [i[0] for i in zip(X.columns,
                                  self.population[g][man][0]) if i[1] == 1]
                score = self.fitnessScore(X[chosen_columns],y,estimator,search_estimator)
                self.population[g][man][1] = score
            passed = time.time()-t
            print "Generation %i took %.2f seconds to fit"%(g,passed) 
            self.nextGeneration(self.population)
        print "Final population is created after %s generations" %n
        
    def mutate(self,men,rate=0.05,bestfitmutants=False):
        n = int(self.genelength*rate)
        g = self.generation
        if bestfitmutants:
            nmen = self.nbestmutants
            for man in np.random.choice(men,size=nmen,replace=False):
                self.population[g][self.manid] = self.population[g-1][man]
                mut = np.random.randint(0,self.genelength-1,n)
                for i in mut:
                    self.population[g][self.manid][0][i] = abs(self.population[g][self.manid][0][i]-1) 
                self.manid += 1
        else:
            nmen = self.nbestcrossmutant
            for man in np.random.choice(men,size=nmen,replace=False):
                mut = np.random.randint(0,self.genelength-1,n)
                self.population[g][self.manid] = self.population[g][man]
                for i in mut:
                    self.population[g][self.manid][0][i] = abs(self.population[g][self.manid][0][i]-1)
                self.manid +=1
    
    def crossover(self,men,returnids=False):
        men2 = list(np.roll(men,1))
        half = self.genelength/2    
        g = self.generation
        for man, man2 in zip(men,men2):
            self.population[g][self.manid]=[self.population[g-1][man][0][0:half],0]
            self.population[g][self.manid][0].extend(self.population[g-1][man2][0][half:half*2])
            self.manid+=1
        if returnids:
            return range(self.manid-len(men),self.manid)

def main():
    
    from sklearn.linear_model import LinearRegression    
    
    df = pd.read_csv("data.csv")
    y = df["y"]
    X = df[set(df.columns)-set(["y"])]
    
    m_lr = LinearRegression()
    evo = EvolutionaryFeatureSubsetting(genelength=1150, populationsize=100)
    evo.fit(X,y,estimator=m_lr,ngenerations=100)
    
    
if __name__ == "__main__":
    main()    
      
