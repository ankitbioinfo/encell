
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def read_data(niche_pred_outdir,Radius):
    confusion_cutoff=0 # no of cell types wants to print
    filename=niche_pred_outdir+'TopCoeff_R'+str(Radius)

    f=open(niche_pred_outdir+'used_CT.txt')
    nameOfCellType={}
    for line in f:
        l=line[0:-1].split('\t')
        nameOfCellType[int(l[0])]=l[1]

    fout=niche_pred_outdir+'niche_prediction_linear/classifier_matrices_'+str(Radius)+'.npz'

    data=np.load(fout,allow_pickle=True)
    coef=data['coef']
    cmn=data['cmn']
    cmn_std=data['cmn_std']
    coef_std=data['coef_std']
    CTFeatures=data['CTFeatures']

    a=np.diag(cmn)
    b=np.diag(cmn_std)
    goodPredictedCellType=np.argsort(-a)
    coeff_cutoff=52

    data={}

    for k in range(len(a)):
        if a[goodPredictedCellType[k]]>=confusion_cutoff:
            meanCoefficients=coef[goodPredictedCellType[k]]
            stdCoefficients=coef_std[goodPredictedCellType[k]]
            highestIndex=np.argsort(-abs(meanCoefficients))

            n=min(coeff_cutoff,len(highestIndex))
            coeff_of_CT=[]
            name_of_the_coeff=[]
            std_of_coeff=[]


            for i in range(n):
            #for i in range(len(highestIndex)):
                l=CTFeatures[highestIndex[i]].split()
                temp=''
                for j in range(len(l)):
                    temp+=nameOfCellType[int(l[j][1:])]
                    if j!=(len(l)-1):
                        temp+='--'
                #print(temp,highestIndex[i],CTFeatures[highestIndex[i]],goodCoefficients[ highestIndex[i]   ])
                integerName=CTFeatures[highestIndex[i]].replace('x','')
                #fw.write(str(highestIndex[i])+'\t'+str('%0.2f'%meanCoefficients[ highestIndex[i]] ) +'\t'+temp+' ('+ integerName  +')\n')
                coeff_of_CT.append(meanCoefficients[ highestIndex[i]])
                name_of_the_coeff.append(temp.replace('_',' '))
                std_of_coeff.append(stdCoefficients[ highestIndex[i]])

            xx=range(len(coeff_of_CT))
            CC=nameOfCellType[goodPredictedCellType[k]]

            data[CC]=[xx, coeff_of_CT, std_of_coeff,name_of_the_coeff]


    return data


def main():
    path=['Healthy/nico_out/','Fibrosis/nico_out/']
    data1=read_data(path[0],0)
    data2=read_data(path[1],0)
    #data3=read_data(path[2],-2)

    mydata=[data1,data2]

    #tname=['juxta(first)','first+second','second']
    tname=['Healthy','Fibrosis']
    CT=sorted(list(data1.keys()))
    topNiche=52

    for i in range(len(CT)):
        allCT=[]
        ally=[]
        for j in range(len(mydata)):
            l=mydata[j][CT[i]]
            #print(j,l)
            mu=l[1]
            std=l[2]
            nc=l[3]
            d={}
            for k in range(len(nc)):
                d[nc[k]]=[mu[k],std[k]]

            nc=l[3][0:topNiche]
            allCT=allCT+nc
            ally.append(d)

        allCT=sorted(list(np.unique(allCT)))
        allCT.remove('Hillock-like')
        #print(allCT)

        ratio={}
        for k in range(len(ally)):
            Y=[]
            Z=[]
            d=ally[k]
            for j in range(len(allCT)):
                value=d[allCT[j]]
                #print(k,value,value[0],value[1])
                Y.append(value[0])
                Z.append(value[1])
            if k==0:
                index=np.argsort(-abs(np.array(Y)))
            Y=np.array(Y)
            ratio[tname[k]]=Y[index]

        allCT=np.array(allCT)
        fig, (ax0) = plt.subplots(1,1, figsize=(10, 2.5))
        df=pd.DataFrame(ratio,index=allCT[index])
        df.plot(ax=ax0,kind='bar',stacked=False,rot=90)
        ax0.set_title(CT[i])
        ax0.legend(loc='upper right',fontsize=7 )  # bbox_to_anchor=(1,0.5)

        savefname=remove_extra_character_from_name(CT[i])

        fig.savefig('./niche_compare/'+savefname+'.png',bbox_inches='tight',transparent=False,dpi=300)


        #print(CT[i],j,allCT)

def remove_extra_character_from_name(name):
    """
    This function removes special characters from the cell type names to avoid throwing an error while saving the figures.
    """
    name=name.replace('/','_')
    name=name.replace(' ','_')
    name=name.replace('"','')
    name=name.replace("'",'')
    name=name.replace(')','')
    name=name.replace('(','')
    name=name.replace('+','p')
    name=name.replace('-','n')
    name=name.replace('.','')
    return name

main()

'''
create_directory(filename)
#for i in range(len(goodPredictedCellType)):
#    print(i,a[goodPredictedCellType[i]])
# top 3 cell type in confusion matrix
for k in range(len(a)):
    if a[goodPredictedCellType[k]]>=confusion_cutoff:
        if  nameOfCellType[goodPredictedCellType[k]] in choose_celltypes:
            flag=1
        else:
            flag=0
        if len(choose_celltypes)==0:
            flag=1

        if flag==1:
            meanCoefficients=coef[goodPredictedCellType[k]]
            stdCoefficients=coef_std[goodPredictedCellType[k]]
            highestIndex=np.argsort(-abs(meanCoefficients))

            n=min(coeff_cutoff,len(highestIndex))
            coeff_of_CT=[]
            name_of_the_coeff=[]
            std_of_coeff=[]

            #fw.write('\n'+str(k+1)+ ' Largest predicted cell type and their top 5 coefficients : '+
            #        nameOfCellType[goodPredictedCellType[k]]+' ( id = '+str(goodPredictedCellType[k])+',  confusion score = '+str('%0.2f'%a[goodPredictedCellType[k]])+')\n')

            for i in range(n):
            #for i in range(len(highestIndex)):
                l=CTFeatures[highestIndex[i]].split()
                temp=''
                for j in range(len(l)):
                    temp+=nameOfCellType[int(l[j][1:])]
                    if j!=(len(l)-1):
                        temp+='--'
                #print(temp,highestIndex[i],CTFeatures[highestIndex[i]],goodCoefficients[ highestIndex[i]   ])
                integerName=CTFeatures[highestIndex[i]].replace('x','')
                #fw.write(str(highestIndex[i])+'\t'+str('%0.2f'%meanCoefficients[ highestIndex[i]] ) +'\t'+temp+' ('+ integerName  +')\n')
                coeff_of_CT.append(meanCoefficients[ highestIndex[i]])
                name_of_the_coeff.append(temp)
                std_of_coeff.append(stdCoefficients[ highestIndex[i]])

            fig,ax=plt.subplots( figsize=figsize)
            xx=range(len(coeff_of_CT))
            yy=np.zeros((len(coeff_of_CT)))
            ax.errorbar(xx, coeff_of_CT, yerr=std_of_coeff,fmt='o',markeredgewidth=0,markerfacecolor=None,markeredgecolor=None,linewidth=1,capsize=1,markersize=2,elinewidth=0.1,capthick=0.75)#markeredgecolor='blue',markerfacecolor='blue',
            ax.plot(xx,yy,'k-',linewidth=0.2)
            #ax.set_ylabel('value of coeff.')
            #ax.set_xlabel('name of the coeff.')
            #titlename=nameOfCellType[goodPredictedCellType[k]]+', conf score = {0:.3f}'.format(a[goodPredictedCellType[k]]) +'$\pm$'+str('%0.3f'%b[goodPredictedCellType[k]])
            titlename=nameOfCellType[goodPredictedCellType[k]]+', conf. score = {0:.3f}'.format(a[goodPredictedCellType[k]]) +'$\pm$'+str('%0.3f'%b[goodPredictedCellType[k]])

            titlename=titlename.replace('_',' ')
            ax.set_title(titlename,fontsize=7)


            ax.set_xticks(xx)
            for ind in range(len(name_of_the_coeff)):
                name_of_the_coeff[ind]=name_of_the_coeff[ind].replace('_',' ')
            ax.set_xticklabels(name_of_the_coeff)
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)
                tick.set_fontsize(7)


            #fig.tight_layout()
            savefname=remove_extra_character_from_name(nameOfCellType[goodPredictedCellType[k]])
            print("The figures are saved: ", filename+'/Rank'+str(k+1)+'_'+savefname+'.'+saveas)
            fig.savefig(filename+'/Rank'+str(k+1)+'_'+savefname+'.'+saveas,bbox_inches='tight',transparent=transparent_mode,dpi=300)
            if showit:
                pass
            else:
                plt.close('all')
'''
