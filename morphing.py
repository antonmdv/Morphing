#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 15:00:43 2017

@author: antonmdv
Author: Anton Medvedev
Project: Morphing
Course: CSE 5280
Prof: Dr Reibiero

"""

from PIL import Image
import numpy
import matplotlib.pyplot as plt
import scipy.misc

im1 = Image.open("P1.jpg")
im2 = Image.open("P2.jpg")
height = 480
width = 340
imSize = [width,height]

numMorphedFrames = 10

#Constants for line weight equation
a = 0.2  
b = 1.25 
m = 0.1  

#11 Lines in total were used

#Src
#P
srcP = list()
srcP.append([200,72])
srcP.append([94,142])
srcP.append([84,142])
srcP.append([306,363])
srcP.append([100,145])
srcP.append([237,190])
srcP.append([131,170])
srcP.append([304,307])
srcP.append([161,137])
srcP.append([204,275])
srcP.append([207,180])
#Q
srcQ = list()
srcQ.append([94,142])
srcQ.append([84,142])
srcQ.append([306,363])
srcQ.append([189,258])
srcQ.append([207,208])
srcQ.append([205,207])
srcQ.append([304,307])
srcQ.append([207,151])
srcQ.append([204,275])
srcQ.append([207,180])
srcQ.append([272,205])

#dest
#P
destP = list()
destP.append([243,55])
destP.append([91,158])
destP.append([65,147])
destP.append([290,393])
destP.append([90,140])
destP.append([250,205])
destP.append([112,172])
destP.append([300,327])
destP.append([161,137])
destP.append([204,275])
destP.append([207,180])
#Q
destQ = list()
destQ.append([91,158])
destQ.append([65,147])
destQ.append([290,393])
destQ.append([200,272])
destQ.append([208,215])
destQ.append([208,212])
destQ.append([300,327])
destQ.append([225,172])
destQ.append([204,275])
destQ.append([207,180])
destQ.append([272,205])


print('Start Morphing')

#returns perpendicular vector
def perpendicular(a) :
    b = numpy.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

#DeltaP and DeltaQ
dP = (numpy.array(destP)-numpy.array(srcP))/(numMorphedFrames+1)
dQ = (numpy.array(destQ)-numpy.array(srcQ))/(numMorphedFrames+1)


for eachFrame in range(0,numMorphedFrames):
    
    #Create new image and get interpolated lines
    morphedIm = Image.new("RGB",imSize,"white")
    interpolatedP = numpy.array(srcP) + dP*(eachFrame+1)
    interpolatedQ = numpy.array(srcQ) + dQ*(eachFrame+1)
    
    #Beier-Neely's multiple-line morphing algorithm
    #for each pixel
    for w in range(0,width):
            for h in range(0,height):
                
                #set counters
                pixel = numpy.array([w,h])
                DSUM1 = numpy.array([0.0,0.0])
                DSUM2 = numpy.array([0.0,0.0])
                weightsum = 0
            
                #for each line
                for line in range(int(numpy.array(srcP).size/2)):
                    
                    #get Points
                    P  = interpolatedP[line,:];
                    P1 = numpy.array(srcP)[line,:];
                    P2 = numpy.array(destP)[line,:];
                    Q  = interpolatedQ[line,:];
                    Q1 = numpy.array(srcQ)[line,:];
                    Q2 = numpy.array(destQ)[line,:];
                    
                    #calculate u and v
                    U =((pixel-P).dot(Q-P))/(numpy.linalg.norm(Q-P)*numpy.linalg.norm(Q-P))
                    V =((pixel-P).dot(perpendicular(Q-P)))/(numpy.linalg.norm(Q-P))
                    
                    #calculate X'
                    xPrime1 = P1 + U*(Q1-P1) + (V*perpendicular(Q1-P1))/numpy.linalg.norm(Q1-P1)
                    xPrime2 = P2 + U*(Q2-P2) + (V*perpendicular(Q2-P2))/numpy.linalg.norm(Q2-P2)
                    
                    #calculate displacement
                    displacement1 = xPrime1 - pixel
                    displacement2 = xPrime2 - pixel
                    
                    #get shortest distance from P to Q
                    if(U >= 1):
                        shortestDist = numpy.linalg.norm(Q-pixel)
                    elif(U <= 0):
                        shortestDist = numpy.linalg.norm(P-pixel)
                    else:
                        shortestDist = abs(V)
                    
                    #calculate line weight
                    lineWeight = ((numpy.linalg.norm(P-Q)**m)/(a+shortestDist))**b
                    
                    #displace sums
                    DSUM1 = DSUM1+lineWeight*displacement1
                    DSUM2 = DSUM2+lineWeight*displacement2
                    weightsum += lineWeight
                
                #displace X' with the sums 
                xPrime1 = pixel +DSUM1/weightsum
                xPrime2 = pixel +DSUM2/weightsum
                
                #get destenation in the new image
                srcX = int(xPrime1[0])
                srcY = int(xPrime1[1])
                destX = int(xPrime2[0])
                destY = int(xPrime2[1])

                #if pixel is in range of the picture,then get color
                #else get color from current pixel
                if(srcX in range(0,width) and srcY in range(0,height)):
                    srcRGB = im1.getpixel((srcX,srcY))
                else:
                    srcRGB = im1.getpixel((w,h))    
                    
                if(destX in range(0,width) and destY in range(0,height)):
                    destRGB = im2.getpixel((destX,destY))
                else:
                    destRGB = im2.getpixel((w,h))
                
                #Cross-disolve weight
                wI2 = float(2*(eachFrame+1)*(1/(numMorphedFrames+1)))
                wI1 = float(2-wI2)
                
                #Cross-dissolve color
                R = (wI1*srcRGB[0] + wI2*destRGB[0])/2
                G = (wI1*srcRGB[1] + wI2*destRGB[1])/2
                B = (wI1*srcRGB[2] + wI2*destRGB[2])/2
                     
                #put color on image
                morphedIm.putpixel((w,h),(int(R),int(G),int(B)))
                
    #save image
    scipy.misc.imsave('morph'+str(eachFrame+1)+'.jpg',morphedIm)
    print('Image '+str(eachFrame+1)+'was saved')
        
            
print('Finish Morphing')