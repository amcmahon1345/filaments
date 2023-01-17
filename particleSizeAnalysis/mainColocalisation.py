import pandas as pd
from Clusterer import Clusterer
import matplotlib.pyplot as plt
import numpy as np

#Parameters
filepath = 'V:\\Virus Group\\Papers\\Vitro Filaments\\Figure 3 SARS\\20201201\\'
redLaserOnFile = 'redPositiveSet3'
greenLaserOnFile = 'greenPositiveSet3'
redChannelForProteinsOfInterest = True
minClusterSize = 50	
dBScanEpsilon = 100
useConvexHullOfPoints = False	
useHierarchicalDBScan = True
maximumSphericalVirionDimension = 300
colocalisationDistance = 200
yFOVSize = 80000

#Read in the data
redLaserData = pd.read_csv(filepath + redLaserOnFile + '.csv')
xyRedLaserData  = redLaserData[['X (nm)','Y (nm)']]
greenLaserData = pd.read_csv(filepath + greenLaserOnFile + '.csv')
xyGreenLaserData  = greenLaserData[['X (nm)','Y (nm)']]

#Use channel info
redXYData = xyRedLaserData[redLaserData['Channel'] == 1]
greenXYData = xyGreenLaserData[greenLaserData['Channel'] == 0]

##Spherical virion analysis
#Cluster channels, colocalise channels, filter virions on colocalisation and size from both channels
if useHierarchicalDBScan == True:
	redClusterLabels = Clusterer.hDBScan(redXYData, minClusterSize)
	greenClusterLabels = Clusterer.hDBScan(greenXYData, minClusterSize) 
	
else:
    redClusterLabels = Clusterer.dBScan(redXYData, minClusterSize, dBScanEpsilon)
    greenClusterLabels = Clusterer.dBScan(greenXYData, minClusterSize, dBScanEpsilon) 
	

redX_centroids,redY_centroids,redMajorAxis,redMinorAxis,redLocalisationsInClusters = Clusterer.fitEllipseAndFindEllipseParametersAndCentroidAndAreaOfClusters(redXYData, redClusterLabels, useConvexHullOfPoints, yFOVSize)
greenX_centroids,greenY_centroids,greenMajorAxis,greenMinorAxis,greenLocalisationsInClusters = Clusterer.fitEllipseAndFindEllipseParametersAndCentroidAndAreaOfClusters(greenXYData, greenClusterLabels, useConvexHullOfPoints, yFOVSize)

colocalisedRedClusterLabels, colocalisedGreenClusterLabels = Clusterer.colocaliseClusters(redX_centroids,redY_centroids,greenX_centroids,greenY_centroids,redClusterLabels,greenClusterLabels,colocalisationDistance)

x_centroids = []
y_centroids = []
majorAxis = []
minorAxis = []
localisationsInClusters = []

if redChannelForProteinsOfInterest == True:
	for i in colocalisedRedClusterLabels:
		for k in i:
			for j in k:
				try:
					j = int(j)
				except:
					print(j)
					print(type(j))
				x_centroids.append(redX_centroids[j])
				y_centroids.append(redY_centroids[j])
				majorAxis.append(redMajorAxis[j])
				minorAxis.append(redMinorAxis[j])
				localisationsInClusters.append(redLocalisationsInClusters[j])
else:
	for i in colocalisedGreenClusterLabels:
		for j in i:
			j = int(j)	
			x_centroids.append(greenX_centroids[j])
			y_centroids.append(greenY_centroids[j])
			majorAxis.append(greenMajorAxis[j])
			minorAxis.append(greenMinorAxis[j])
			localisationsInClusters.append(greenLocalisationsInClusters[j])

results = pd.DataFrame({'x_centroids': x_centroids})
results.insert(1,'y_centroids', y_centroids)
results.insert(2,'Number of localisations in the cluster', localisationsInClusters)
results.insert(3,'Major Axis', majorAxis)
results.insert(4,'Minor Axis', minorAxis)
results.to_csv(redLaserOnFile + 'LocalisationsInClusters.csv')
