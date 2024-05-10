
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#print working directory
import os
print(os.getcwd())

figure_counter=1

#point to where TL output files are kept
FILEBASE='/Users/christopherjakobson/Dropbox/CONDENSeq Background Data/'
#make figure output directory if not already present
if not os.path.exists(FILEBASE + 'analyze_counts_figures'):
    os.makedirs(FILEBASE + 'analyze_counts_figures')

#import a list of text files
files_to_import=['2022_10_31_ReadCountAnalysis_Kmeans8clusters_thresholded_fullTable.txt',
                '2022_11_28_ReadCountAnalysis_Kmeans8clustersB_thresholded_fullTable.txt',
                '2022_11_28_ReadCountAnalysis_Kmeans8clustersC_thresholded_fullTable.txt']

arm_labels=['Gal=>Raf','Raf=>Raf','Raf=>no selection']
day_labels=['day0','day2','day5','day7']

#import the data in files_to_import as a list of dataframes
dataframes = [pd.read_csv(FILEBASE + file, sep='\t') for file in files_to_import]

#check the first few rows of the first dataframe
#print(dataframes[0].head())

#store columns 6 to 9 of each dataframe as a list of numpy arrays
counts = [df.iloc[:, 5:9].values for df in dataframes]

#check the shape of the first array
#print(counts[0].shape)

#remove low counts
count_thresh=1
#set values in counts below 1 to 0
counts = [np.where(counts[i] < count_thresh, 0, counts[i]) for i in range(len(counts))]      

#print(counts[0])

#normalize columns to a total of 1e6 reads
total_reads = [np.sum(counts[i], axis=0) for i in range(len(counts))]
#print(total_reads[0])
counts = [counts[i]/total_reads[i][None,:]*1e6 for i in range(len(counts))]

#print(counts[0])

#take the mean of the three replicates
means = [np.zeros((int(len(counts[i])/3),4)) for i in range(len(counts))]
#print(means[0].shape)

for i in range(len(counts)):
    for j in range(len(day_labels)):
        temp_mat=np.zeros((int(len(counts[i])/3),2)) #initialize a temporary matrix to store the means
        for k in range(0,2):
            temp_mat[:,k]=counts[i][k+0:-1:3,j]
        means[i][:,j]=np.mean(temp_mat, axis=1)

#print(means[0])


#fill in pseudocounts in mean
read_thresh=0
pseudo_count=0.1

for i in range(len(means)):
    #zeros get pseudocount if the gene ever shows up in that arm
    v_ever_observed=np.sum(means[i]>read_thresh, axis=1)>0
    #print(v_ever_observed)
    means[i][means[i]<=read_thresh]=pseudo_count
    #nan missing rows
    means[i][~v_ever_observed,:]=np.nan
    #remove stray zeros
    means[i][means[i]==0]=np.nan

#print(means[0])

#cross-plot
subplot_counter=1
plt.figure()
for i in range(len(means)):
    for j in range(4):
        for k in range(j+1,4):
            plt.subplot(3,6,subplot_counter)
            plt.scatter(means[i][:,j], means[i][:,k], s=0.5)
            plt.xlabel(day_labels[j])
            plt.ylabel(day_labels[k])
            plt.title(arm_labels[i])
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim(1e-2,1e4)
            plt.ylim(1e-2,1e4)
            plt.gca().set_aspect('equal', adjustable='box')
            #add spacing between subplots
            plt.subplots_adjust(wspace=0.75, hspace=0.75)
            subplot_counter=subplot_counter+1
#make fullscreen
mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
#save figure as svg
#plt.show()
plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.svg', format='svg', dpi=300)
plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.png', format='png', dpi=300)
figure_counter=figure_counter+1
#plt.close()


#renormalize to day 0 for each arm
means_norm = [np.zeros((int(len(counts[i])/3),4)) for i in range(len(counts))]

for i in range(len(means)):
    for j in range(4):
        means_norm[i][:,j]=means[i][:,j]/means[i][:,0]

#print(means_norm[0])


#take sum of log10 of each row to use as "score"
scores = [np.zeros((int(len(counts[i])/3),1)) for i in range(len(counts))]

for i in range(len(means)):
    scores[i]=np.sum(np.log10(means_norm[i]), axis=1)

#print(scores[0])

#plot histograms of scores
plt.figure()
for i in range(len(scores)):
    plt.subplot(2,3,i+1)
    plt.hist(scores[i], bins=100, alpha=0.5, label=arm_labels[i])
    plt.title(arm_labels[i])
    plt.xlabel('score')
    plt.ylabel('frequency')
mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
#save figure as svg
#plt.show()
plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.svg', format='svg', dpi=300)
plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.png', format='png', dpi=300)
figure_counter=figure_counter+1
#plt.close()


#sort scores and plot rows as a heatmap in ascending order
for i in range(len(scores)):
    sorted_indices = np.argsort(scores[i])[::-1]
    plt.figure()
    plt.imshow(means_norm[i][sorted_indices,:], aspect='auto',vmin=0, vmax=10, cmap='seismic')
    plt.colorbar()
    plt.title(arm_labels[i])
    plt.xlabel('day')
    plt.ylabel('gene')
    plt.yticks([])
    plt.xticks([0,1,2,3], day_labels)
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    #save figure as svg
    #plt.show()
    plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.svg', format='svg', dpi=300)
    plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.png', format='png', dpi=300)
    figure_counter=figure_counter+1
    #plt.close()

    #plot some example time courses of means_norm on the same axes
    plt.figure()
    for j in range(0,len(sorted_indices),1000):
        plt.plot(means_norm[i][sorted_indices[j],:])
    plt.title(arm_labels[i])
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    #save figure as svg
    #plt.show()
    plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.svg', format='svg', dpi=300)
    plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.png', format='png', dpi=300)
    figure_counter=figure_counter+1
    #plt.close()




#match to vidal ids, uniprot ids, and biophysical properties
#vidal ids are the same for all three dataframes
v_temp = dataframes[0]["VidalID"]
vidal_ids=v_temp[0:-1:3].values


#import vidal dictionary from /Users/christopherjakobson/Dropbox/CONDENSeq Background Data/CONDENSeqBinAnalysis.xlsx
#Arm B tab has all the vidal IDs
vidal_dict = pd.read_excel(FILEBASE + 'CONDENSeqBinAnalysis.xlsx', sheet_name='Arm B', index_col=0)

#print(vidal_dict.head())

#get rows of vidal_dict for which Clone matches vidal_ids
gene_info=vidal_dict[vidal_dict.index.isin(vidal_ids)]

print(gene_info.head())


#scatter the scores against one another
plt.figure()
subplot_counter=1
for i in range(len(scores)):
    for j in range(i+1,len(scores)):
        plt.subplot(2,3,subplot_counter)
        plt.scatter(scores[j], scores[i], s=0.5)
        plt.xlabel(arm_labels[j])
        plt.ylabel(arm_labels[i])
        plt.title('scores')
        plt.xlim(-10,10)
        plt.ylim(-10,10)
        plt.gca().set_aspect('equal', adjustable='box')
        #add spacing between subplots
        plt.subplots_adjust(wspace=0.75, hspace=0.75)
        subplot_counter=subplot_counter+1

mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
#save figure as svg
#plt.show()
plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.svg', format='svg', dpi=300)
plt.savefig(FILEBASE + 'analyze_counts_figures/analyze_counts_figure' + str(figure_counter) + '.png', format='png', dpi=300)
figure_counter=figure_counter+1
#plt.close()


#categorize based on scores
templating_nontoxic_index=np.where((scores[1]>0)*(scores[1]<=scores[0]))[0]
templating_toxic_index=np.where((scores[1]>0)*(scores[1]>scores[0]))[0]
aggregating_index=np.where((scores[1]<=0)*(scores[0]>0.5))[0]
#all_other_index=np.where((scores[1]<=0)*(scores[0]<=0.5))[0]


#print(templating_nontoxic_index)

#assign categories
v_categories=list(range(len(scores[0])))
for i in range(len(v_categories)):
    v_categories[i]='all_other'
for i in templating_nontoxic_index:
    v_categories[i]='templating_nontoxic'
for i in templating_toxic_index:
    v_categories[i]='templating_toxic'
for i in aggregating_index:
    v_categories[i]='aggregating'

print(v_categories[0:9])


biophys_dict=pd.read_csv(FILEBASE + 'BiophysicalFeatures_2022_05_23.txt', sep='\t')

biophys_info=biophys_dict[biophys_dict["Vidalno"].isin(vidal_ids)]

print(biophys_info.head())




