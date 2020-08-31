#!/home/dmollet/.conda/envs/python_env/bin/python3

#functions to be imported !

import glob
import os
import sys
import csv
import collections
from PIL import Image
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler, Normalizer, normalize
import operator
import numpy as np
from Bio import SeqIO
from dict_of_codons import * 
import scipy
import pygmnormalize as gmn
import cv2
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

def get_codon_usage(data, file_list) :

    trans = str.maketrans('ATGC', 'TACG')

    all_images = []

    #for each strain
    for num,i in enumerate(data) :

        str_seq = str(i)
        
        #get stats file
        f = "stats_all"+file_list[num][5:-5]+"stats_all"
        
        print(f)
        gene_positions = []
        number_of_genes = 0

        #put it into a list
        with open(f) as s:
            stat_file = csv.reader(s, delimiter=' ')

            #number_of_genes = sum(1 for row in stat_file) - 1

            next(stat_file, None)         #header

            for j in stat_file :
                if(j[1] != '0') :
                    gene_positions.append([int(j[2]), int(j[3]), int(j[4]), int(j[1])])

        #sort
        gene_positions.sort(key=lambda x: int(x[3]))

        number_genes = len(gene_positions)

        #648 entered by hand ! find a better way
        image = np.zeros(shape=(648, 128))

        #for each gene
        for pos,k in enumerate(gene_positions) :

            codons = []

            direction = k[2]*2 - 1 #1 if positive, -1 if negative

            #get all the codons
            if k[2] == 1 :
                codons = [str_seq[j:j+3] for j in range(k[0]-1, k[1]-2, 3)]

            elif k[2] == 0:
                codons = [str_seq[j-3:j][::-1].translate(trans) for j in range(k[1], k[0]+1, -3)]

            #count em
            temp_dic = CodonsDictWide.copy()

            for codon in codons :
                c = codon.upper()
                try: 
                    if direction == 1 :
                        temp_dic[c] += 1 
                    elif direction == -1 :
                        temp_dic['n'+c] += 1
                except:
                    pass

            #in an image, order is same as initial dict
            image[pos,:] = np.fromiter([temp_dic[i] for i in sorted(temp_dic.keys())], dtype=float)
            #print(direction)
            #print(image[pos,:])

        #---------------- TESTING
        if False:

            #minmax scaling
            scaler = MinMaxScaler(feature_range = (0,255))
            minmax_image = scaler.fit_transform(image[:number_genes,...]).astype(np.uint8)

            #----------------
            #tmm
            tmm_image = gmn.tmm_normalization(minmax_image)

            #----------------
            #histogram equalization
            equ_image = cv2.equalizeHist(minmax_image)

            #----------------
            #contrast limited adaptive hist equ
            clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
            clahe_image = clahe.apply(minmax_image)
            #clahe_image = clahe.apply(image[:number_genes,...].astype(np.uint8))

            #----------------

            [minmax_image, tmm_image, equ_image, clahe_image] = [np.concatenate((i, image[number_genes:,...].astype(np.uint8)), axis=0) for i in [minmax_image, tmm_image, equ_image, clahe_image]]

            sns.set_style('darkgrid')

            all_titles = ['raw', 'minmax', 'tmm', 'equ hist', 'clahe']
            all_images = [image, minmax_image, tmm_image, equ_image, clahe_image]
            all_colors = ['black', 'red', 'green', 'blue', 'orange']

            #----------------
            #combine all images
            for i in all_images :
                i = Image.fromarray(i).convert('L')

            fig, axes = plt.subplots(nrows=1, ncols=5)
            fig.tight_layout()

            for i in range(0,len(all_images)):
                ax = fig.add_subplot(1,len(all_images)+1,i+1)
                plt.imshow(all_images[i])#, cmap="Greys")
                plt.axis('off')
                axes[i].axis('off')
                ax.set_title(all_titles[i], fontsize=8)

            #fig.savefig('all_images.jpg', bbox_inches='tight')
            fig.savefig('images/widecomp_'+file_list[num][6:-5]+'jpg', bbox_inches='tight')

            #hists
            #fig = plt.figure(figsize=(20,20))
            fig, ax = plt.subplots(figsize=(20,20))

            for i in range(0,len(all_images)):

                #f = sns.distplot(all_images[i][:number_genes,...].flatten(), color=all_colors[i], kde=False, fit=stats.skewnorm, fit_kws={'color':all_colors[i]}, label=all_titles[i], hist=None)
                f = sns.distplot(all_images[i][:number_genes,...].flatten(), color=all_colors[i], bins=200, kde=False, fit=stats.skewnorm, fit_kws={'color':all_colors[i]}, label=all_titles[i])
                f.set_xlim(-10,255)
                f.set_ylim(0,0.15)

                l = f.lines[i]
                x = l.get_xydata()[:,0]
                y = l.get_xydata()[:,1]
                f.fill_between(x,y, color=all_colors[i], alpha=0.3)

            plt.axis('on')
            #plt.axvline(x=127, ymin=0, ymax=0.035, color='red')
            ax.legend(ncol=2, loc='center left', fontsize=15, title='Scaling methods', title_fontsize=20)
            #fig.savefig('all_hists.jpg', bbox_inches='tight')
            fig.savefig('images/widehist_'+file_list[num][6:-5]+'jpg', bbox_inches='tight')

            #break

        #----------------
        #final image

        clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
        clahe_image = clahe.apply(image[:number_genes,...].astype(np.uint8))

        final_image = np.concatenate((clahe_image, image[number_genes:,...].astype(np.uint8)), axis=0)

        print(final_image.shape)

        #final_image = Image.fromarray(final_image).convert('L')
        #final_image.save('images/wide_'+file_list[num][6:-5]+'png')
        #print("done !")

        all_images.append(final_image)

        #break

    return(all_images)
