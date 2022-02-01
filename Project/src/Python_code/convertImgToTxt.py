import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#import otsuthreshold as ot
import scipy.misc


def deleteIllegalOnes(img):
    """
    Funktion gets an Array in geometry format and deletes all ones (obstaclecells),
    that contain fluidcells on opposite sites
    """
    corrections=1
    while (corrections>0):
        corrections=0
        for i in range(1,img.shape[0]-1):
            for j in range(1,img.shape[1]-1):
                if img[i,j]==1:
                    neighbourisfluid=[img[i,j+1]==0,img[i+1,j]==0,img[i,j-1]==0,img[i-1,j]==0]
                    if ((neighbourisfluid[0] and neighbourisfluid[2]) or (neighbourisfluid[1] and neighbourisfluid[3]))  :
                        img[i,j]=0
                        corrections=corrections+1
        for i in range(1,img.shape[0]-1):
            for j in [0,img.shape[1]-1]:
                if img[i,j]==1:
                    if (img[i+1,j]==0 and img[i-1,j]==0):
                         img[i,j]=0
                         corrections=corrections+1
        for j in range(1,img.shape[1]-1):
            for i in [0,img.shape[0]-1]:
                if img[i,j]==1:
                    if (img[i,j-1]==0 and img[i,j+1]==0):
                        img[i,j]=0
                        corrections=corrections+1
                   
def printToFile(titel,data):
    """
    Gets an string and an numpy Array of Datas an prints the dataas integers 
    into an file with the name contained in the string
    the format ist adapted to the needs of the c++ program
    """
    open(titel+".txt", 'w').close()
    f = open(titel+".txt", "a")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            f.write("{:n},".format(data[i,j])) 
        f.write("\n")     
    f.close()
     
def convertToBinaryGeometry(data):
    """
    Receives an image as numpy array and transforms it to 
    a binary file with 0 for white pixels and 1 for all others
    """
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i,j]==255:
                data[i,j]=0
            else:
                data[i,j]=1    
            
def doublePixelHorizontal(img):
    """
    Doubles the pixelnumber in the horizontal direction and 
    adapts the data to the larger format
    """
    newimage=np.zeros((img.shape[0],2*img.shape[1]))    
    for i in range(img.shape[0]):                  
        for j in range(img.shape[1]):
            newimage[i,2*j]=img[i,j]
            newimage[i,2*j+1]=img[i,j]
    return newimage        

def doublePixelVertical(img):
    """
    Doubles the pixelnumber in the vertical direction and 
    adapts the data to the larger format
    """
    newimage=np.zeros((2*img.shape[0],img.shape[1]))    
    for i in range(img.shape[0]):                  
        for j in range(img.shape[1]):
            newimage[2*i,j]=img[i,j]
            newimage[2*i+1,j]=img[i,j]
    return newimage        

def quadrupelPixels(img):
    """
    Doubles the pixelnumber inboth directions and 
    adapts the data to the larger format
    """
    doubleimg=doublePixelHorizontal(img)
    quadrupelimg=doublePixelVertical(doubleimg)
    return quadrupelimg
    


#give the path to the image
imgpath = "test.jpg"

#converts the image into a 3 dimensinal array (may be important for inflow outflow, different bC etc)
img = np.array(mpimg.imread(imgpath)) 
img.astype(np.uint8)

#converts the data to a grey scale i.e. black white image
blackwhitedata = np.mean(img,axis=2).astype(np.uint8)

# creates a binary array with only zeros and ones for fluid and obstacles respectively
convertToBinaryGeometry(blackwhitedata)

#deletes all obstacles that contain fluidcells on opposite sites
deleteIllegalOnes(blackwhitedata)

#prints the data to a txt file
printToFile(imgpath,blackwhitedata)



#blackwhitedata= quadrupelPixels(blackwhitedata)
     
  
