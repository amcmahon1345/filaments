import pandas as pd
from Clusterer import Clusterer
import matplotlib.pyplot as plt
import numpy as np
import os

#Parameters

#RESULTS ARE CURRENTLY NOT SAVED- LINE IS COMMENTED OUT AT BOTTOM

filepath = 'V:\\Virus Group\\Papers\\Vitro Filaments\\Figure 3 Spheres\\Figure 3 SARS\\2020815\\'
minClusterSize = 200
dBScanEpsilon = 30
useConvexHullOfPoints = False	
useHierarchicalDBScan = False
useRedChannel = False	
useGreenChannel = True
maximumSphericalVirionDimension = 300
colocalisationDistance = 50
yFOVSize = 80000

#Read in the data
	
fileArray = []
for d,r,f in os.walk(filepath):
	for file in f:
		if file.endswith(".csv"):
			fileArray.append(os.path.join(d,file))
			
for fullFilePath in fileArray:
			laserData = pd.read_csv(fullFilePath)
			xyLaserData  = laserData[['X (nm)','Y (nm)']]

			#Segment into channels
			if useRedChannel == True:
				redXYData = xyLaserData[laserData['Channel'] == 1]
			if useGreenChannel == True:
				greenXYData = xyLaserData[laserData['Channel'] == 0]

			##Spherical virion analysis
			#Cluster channels, colocalise channels, filter virions on colocalisation and size from both channels
			if useHierarchicalDBScan == True:
				if useRedChannel == True:
					redClusterLabels = Clusterer.hDBScan(redXYData, minClusterSize)
				if useGreenChannel == True:
					greenClusterLabels = Clusterer.hDBScan(greenXYData, minClusterSize) 
				
			else:
				if useRedChannel == True:
					redClusterLabels = Clusterer.dBScan(redXYData, minClusterSize, dBScanEpsilon)
				if useGreenChannel == True:
					greenClusterLabels = Clusterer.dBScan(greenXYData, minClusterSize, dBScanEpsilon) 
					
			if useRedChannel == True:	
				#Clusterer.plotClusters(redXYData, redClusterLabels, yFOVSize)
				x_centroids,y_centroids,majorAxis,minorAxis,localisationsInClusters = Clusterer.fitEllipseAndFindEllipseParametersAndCentroidAndAreaOfClusters(redXYData, redClusterLabels, useConvexHullOfPoints, yFOVSize)
			if useGreenChannel == True:
				#Clusterer.plotClusters(greenXYData, greenClusterLabels, yFOVSize)
				x_centroids,y_centroids,majorAxis,minorAxis,localisationsInClusters = Clusterer.fitEllipseAndFindEllipseParametersAndCentroidAndAreaOfClusters(greenXYData, greenClusterLabels, useConvexHullOfPoints, yFOVSize)

			results = pd.DataFrame({'x_centroids': x_centroids})
			results.insert(1,'y_centroids', y_centroids)
			results.insert(2,'Number of localisations in the cluster', localisationsInClusters)
			results.insert(3,'Major Axis', majorAxis)
			results.insert(4,'Minor Axis', minorAxis)
			
			head_tail = os.path.split(fullFilePath)
			
			if useRedChannel == True:
				results.to_csv(head_tail[1] + 'Red' + 'LocalisationsInClusters.csv')
			if useGreenChannel == True:
				results.to_csv(head_tail[1] + 'Green' + 'LocalisationsInClusters.csv')
