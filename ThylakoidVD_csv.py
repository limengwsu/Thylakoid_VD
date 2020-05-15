# -*- coding: utf-8 -*-
"""
ThylakoidVD_csv.py was written to analyze the vertical dimension of thylakoid grana
stacks, with input of X,Y txt/csv files, output a file of calculated Repeat distances,
Stroma gap width... and the plots of data and fit to inspect the quality of fittings.
Created by Li, Meng @WSU in 2017, modified a few times
email: m.li (at) wsu.edu
This was the first meaningful code I wrote for photosynthesis research, which
has many newbie charateristics, but it worked for our purpose.

NOTE: this code may not work if the lumen is resolved and darkly stained
"""
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages#save multi fig to pdf

def func(x, u, d):
    return 1-np.exp((-(x-u)**2)/(2*d**2))# define desired curve function after processing
#function to return how many columns of data of a file
def filewidth(datafile):#can be used after open a file
    line=datafile.readline()
    while len(line.split())==0:
        line=datafile.readline()
    linelen = len(line.split())
    return linelen

#The following function finds a segment of peak to peak, normalizes and fits it to the given function.
def PtoPprocess(x,y,i,nm):
    #input of nm is the number of data points = 1 nm, which should be a float
    #i is expected to be 3, where initial peak is guessed, pass from ixy
    xnm = int(6*nm) #a peak needs be the maximum within the following 6 nm (4 nm membrane + 2 nm lumen)
    min_half_Sgap = 3.0 #the min of (stroma_gap+ 1mbrane)/2
    while y[i]< max(y[i-3:i+xnm]):#x times nm defines the peak within x nm
        i = i+1
    #ymax = y[i] #ymax can be used instead of yL & YR if stroma gap is not of interest
    j = i
    while y[j]>min(y[j-3:j+xnm]):
        j = j+1
    ymin = y[j]
    Lshder = j-1# try to find left shoulder
    while y[Lshder] > y[Lshder+1] or x[j]-x[Lshder]<min_half_Sgap:
        #assuming stroma gap at least 2 nm
        Lshder = Lshder - 1
    yL = y[Lshder+1]#Left shoulder value
    uguess = x[j]
    yPtoV = np.array(y[i:j])#creat an array of y from peak(P)to valley(V), here peak is shoulder peak
    ynorm = (yPtoV-ymin)/(yL-ymin)#shder as 1 normalize
    xPtoV = x[i:j]#creat a list of x from peak(P)to valley(V)

    #from valley to peak---find peaks and valleys
    k = j
    while y[k]< max(y[k-3:k+xnm]):
        k = k+1
    #ymax = y[k]
    Rshder = j+1#Right shoulder x index (now scan from j to k)
    while y[Rshder] > y[Rshder-1] or x[Rshder]-x[j]<min_half_Sgap:
        #assuming P-to-P at least 6 nm
        Rshder = Rshder + 1
    yR = y[Rshder-1]#Right shoulder value
    dguess = (x[Rshder-1]-x[Lshder+1])/2#using shoulder width to estimate d
    yVtoP = np.array(y[j:k+1])#creat an array of y from Valley(V) to Peak(P)
    ynorm2 = (yVtoP-ymin)/(yR-ymin)
    xVtoP = x[j:k+1]#Creat a list of x from Valley(V) to Peak(P)
    i = k
    xnxs = np.concatenate((xPtoV, xVtoP)) # generate a peak-to-peak(2nd P omitted) xnormalized axis
    ynxs = np.concatenate((ynorm, ynorm2))#np.array(ynorm.tolist()+ynorm2.tolist())
    plt.plot(xnxs, ynxs, 'bo')

    # Fit for the parameters u,d of the function `func` using shoulder to shoulder data
    xshder = np.array(x[Lshder+1:Rshder])
    yLS = np.array(y[Lshder+1:j])
    ynormLS = (yLS-ymin)/(yL-ymin)
    yRS = np.array(y[j:Rshder])
    ynormRS = (yRS-ymin)/(yR-ymin)
    yshder = np.concatenate((ynormLS, ynormRS))#np.array(ynormLS.tolist()+ynormRS.tolist())
    popt, pcov = curve_fit(func, xshder, yshder,bounds = (0,[2*uguess,2*dguess]))
    plt.plot(xshder, func(xshder, *popt), 'r-')

    return(popt,i,xnxs,ynxs)

#from x,y data, nomarlize, fit and plot, NOTE that PtoPprocess is call here as well.
def getud(xdata,ydata,width=4):
    #process data and plot; width is the number of data points within ~1 nm
    ixy = 3
    uf = [] # a list of u, peak postions, f for function, distinguish later u list
    df = [] # a list of d, the sigma, FWHM/2.355
    #iterate the function of finding Peaks and Valleys,
    #plot as well, return valley postions and sigma
    global udcalc
    try:
        while ixy>1:
            udcalc = PtoPprocess(xdata,ydata,ixy,width) #udcalc is the returned tuple
            ixy = udcalc[1] # update i
            udarray = udcalc[0] #optimized u,d in array format
            uf = uf + [udarray.tolist()[0]]
            df = df + [udarray.tolist()[1]]
    except IndexError:
        #replot last set fit data to show label
        plt.plot(udcalc[2],udcalc[3], 'bo', label='data')
        #replot last set fit to show label
        plt.plot(udcalc[2], func(udcalc[2], *udcalc[0]), 'r-', label='fit')
    except RuntimeError:
        print('Cannot fit certain data, RuntimeError')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    return(uf,df)

#using u and d to calculate the lumen/stromal gap widths and repeat distance
def CalcLSRd(u,d):
    #calculate the repeat distance and Stromal gap distances
    rpd = []
    lumen = []
    ncalc = len(u)#the length of u equals how many times of calculation P-V-P
    print(ncalc)
    if ncalc in [0,1]:
        Sgap = np.array([])
    elif (u[-1]-u[0])/(ncalc-1)<12:#judge if the lumen is resolved
        #the following assumes the first is lumen,
        #the gap to gap is used for repeat distance
        idx = 1
        stromalgap = []
        while idx+2 < ncalc:
            rd = u[idx+2]-u[idx]#calculate V to V distance
            rpd = rpd +[rd]
            lumen = lumen + [d[idx-1]*2.355]#lumen width FWHM
            stromalgap = stromalgap + [d[idx]*2.355]#stromal gap FWHM
            idx = idx +2
        if idx < ncalc:
            lumen = lumen + [d[idx-1]*2.355]#lumen width FWHM
            stromalgap = stromalgap + [d[idx]*2.355]#stromal gap FWHM
        Sgap = np.array(stromalgap)#make this list into array
    else:#this is the situation that only gap is resovled
        Sgap = np.array(d)*2.355#stromal gap array
        idx = 0
        while idx+1 < ncalc:
            rd = u[idx+1] - u[idx]#calculate V to V distance
            rpd = rpd + [rd]
            idx = idx + 1
    lumen = np.array(lumen)
    return(lumen,Sgap,rpd)

#file_name = '/Users/sisixiong/Desktop/HPF_Analysis2020/201709_L_xy.csv'
def process_a_file(file_name):
    #with input of a filename, output files of plots of data and fit,
    #and calculated thylakoid/grana vertical dimension data
    separator = ''
    with open(file_name, 'r') as f:
        contents = f.read().splitlines()
        if '\t' in contents[0]:
            separator = '\t'
        elif ',' in contents[0]:
            separator =','
        else:
            print("CALCULATION STOPPED!!!\n\
                  ERROR IS INEVITABLE IF INPUT FILE IS INCORRECT.\n\
                  PLEASE MAKE SURE THE INPUT FILE HAS CORRECT FORMAT WITH TAB OR COMMA AS SEPARATOR\DELIMITER")
    if separator != '':
        for i, aline in enumerate(contents):
            contents[i] = re.split(separator, aline)
            #use the separator depending on the how the file column is separated
        nsets = int(len(contents[0])/2)#how many sets/measurements of data
        print(nsets, 'sets of data will be anaylized for ' + file_name)
        a = np.array(contents[1:])# remove first row and asign data to np array
        a[a==''] = np.nan#fill empty spaces
        a = a.transpose().astype(float)

        pp = PdfPages(file_name.replace('.csv',' plots.pdf'))
        LSRfile = open(file_name.replace('.csv',' LSRd.txt'),'w')
        for theset in range(0,nsets*2,2):
            #the following block extracts the x,y data at current column/set
            nans = np.where(np.isnan(a[theset]))#giving a tuple of an array (array([...]),)
            if nans[0].shape == (0,):#which means the original column does not have empty space
                the_end = -1
            else:
                the_end = nans[0][0]
            xflt = a[theset][:the_end]
            yflt = a[theset+1][:the_end]
            #the above block extracts the x,y data at current column/set

            one_nm = 1/(xflt[1]-xflt[0]) #calculates how many x points needed for one nanometer
            #the following plot the original data
            plt.title('data'+str(theset/2+1))
            plt.plot(xflt,yflt,'g-')
            pp.savefig()
            plt.close()
            #The following process one set of x, y data,
            #calculating the interested values and save the fit picture
            plt.figure(theset/2+1)#each time a different figure number
            plt.title('fit '+str(theset/2+1))
            udncalc=getud(xflt,yflt,one_nm)#this is not simple,
            #it also calls the PtoP process function, with plot done at the same time.
            pp.savefig()#save current figure
            plt.close()#close the plot otherwise multiplot on one figure
            LSRinfo = CalcLSRd(*udncalc)#Calculate the Lumen, Stromal Gap, Repeat Distance

            #The following segment make the data to list and then strings for txt file writing.
            uc = udncalc[0]#call the u (calculated) list
            dc = udncalc[1]#call the d (calculated) list
            #ncalc=udncalc[2]# call the ncalc, number of times of calculation
            #print(uc,dc,*LSRinfo)#LSRinfo has lumen, Stgap, Repeat distance in array, array, list format
            ucstr = '\t'.join(str(e) for e in uc)
            dcstr = '\t'.join(str(f) for f in dc)
            lm = LSRinfo[0].tolist()#get lumen list from array
            Stgap = LSRinfo[1].tolist()#Stromal gap
            Rd = LSRinfo[2]#assign repeat distance to Rd
            lmstr = '\t'.join(str(l) for l in lm)
            Ststr = '\t'.join(str(s) for s in Stgap)
            Rdstr = '\t'.join(str(r) for r in Rd)
            #The following write the txt which can be easily transferred to excel sheet.
            LSRfile.write('Valley Centers '+str(theset/2+1)+'\t'+ucstr+'\n')
            LSRfile.write('sigmas '+str(theset/2+1)+'\t'+dcstr+'\n')
            LSRfile.write('lumen width '+str(theset/2+1)+'\t'+lmstr+'\n')
            LSRfile.write('Stromal gap width '+str(theset/2+1)+'\t'+Ststr+'\n')
            LSRfile.write('Repeat distances '+str(theset/2+1)+'\t'+Rdstr+'\n')
        #after processing each file, close the file
        pp.close()
        LSRfile.close()
#process_a_file('/Users/sisixiong/Desktop/HPF_Analysis2020/201709_L_xy.csv')

list_of_files = glob.glob('./*.csv')# create the list of file
print(list_of_files)
for afile in list_of_files:
    print('Now processing ', afile)
    process_a_file(afile)
