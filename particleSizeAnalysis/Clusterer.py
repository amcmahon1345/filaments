import numpy as np
import matplotlib.pyplot as plt
import hdbscan
import pandas as pd
import seaborn as sns
import scipy
from itertools import combinations
import random
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import confidence_ellipse

class Clusterer:
    
	def square_distance(x,y): return sum([(xi-yi)**2 for xi, yi in zip(x,y)])
	
	#Shortened version of confidence ellipse just finding the size of the ellipse
	def returnHeightAndWidthOfEllipse(xCoordinates, yCoordinates, n_std):
		cov = np.cov(xCoordinates, yCoordinates)
		pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
		# Using a special case to obtain the eigenvalues of this
		# two-dimensionl dataset.
		ell_radius_x = np.sqrt(1 + pearson)
		ell_radius_y = np.sqrt(1 - pearson)
		ellipse_x = ell_radius_x * np.sqrt(2) * n_std * np.sqrt(cov[0, 0] + cov[1, 1])	#sqrt(2) takes into account rotation needed and radius to diameter
		ellipse_y = ell_radius_y * np.sqrt(2) * n_std * np.sqrt(cov[0, 0] + cov[1, 1])
		
		return ellipse_x,ellipse_y
		
    
	def fitEllipseAndFindEllipseParametersAndCentroidAndAreaOfClusters(X, clusterLabels, useConvexHullOfPoints, yFOVSize):
		# Find the area and centroid and major and minor axis distances. If useConvexHullOfPoints, the ellipse is fit to the convex hull. Otherwise it is fit to all points. The fitting is sqrt of principle moments from the focus.
		
		uniqueLabels = set(clusterLabels)
		
		x_centroids = []
		y_centroids = []
		majorAxis = []
		minorAxis = []
		localisationsInClusters = []
		
		#Go through all of the clusters
		for label in uniqueLabels:
			if label!= -1 & len(X[clusterLabels == label])>2:
				if useConvexHullOfPoints == True:
					hull = scipy.spatial.ConvexHull(X[clusterLabels == label])
					vertices = hull.vertices
					vertexCoordinates = hull.points[vertices, :]
					xCoordinates = np.array([p[0] for p in vertexCoordinates])
					yCoordinates = np.array([p[1] for p in vertexCoordinates])
					try:
						ellipse = EllipseModel()
						ellipse.estimate([[xCoordinates], [yCoordinates]])
						xc, yc, height, width, theta = ellipse.para
					except:
						pass
				else:
					xCoordinates = np.array(X[clusterLabels == label]['X (nm)'])
					yCoordinates = np.array(X[clusterLabels == label]['Y (nm)'])
					
					#These lines fit best ellipses using the original program and plot them
					plt.plot(xCoordinates, yCoordinates,'o', markerfacecolor='g', markeredgecolor='k', markersize=5)
					ellipse = confidence_ellipse.confidence_ellipse(xCoordinates, yCoordinates, plt.gca(), n_std=2.0, edgecolor='red')
					#plt.show()
					
					height,width = Clusterer.returnHeightAndWidthOfEllipse(xCoordinates,yCoordinates, n_std=2.0)
					
					#Remove ones where fit has incorrectly found ellipse
					if isinstance(height, complex) or isinstance(width, complex) or height<0 or width<0: 		
						continue		
					
					else:
						x_centroids.append(np.mean(X[clusterLabels == label].iloc[:,0]))
						y_centroids.append(yFOVSize - np.mean(X[clusterLabels == label].iloc[:,1]))
						localisationsInClusters.append(len(X[clusterLabels == label]))
						
						
						if height > width:
							majorAxis.append(height)
							minorAxis.append(width)
						else:
							majorAxis.append(width)
							minorAxis.append(height)
		return 	x_centroids,y_centroids,majorAxis,minorAxis,localisationsInClusters
		
	def filterClusters():
		return
		
	def colocaliseClusters(redX_centroids,redY_centroids,greenX_centroids,greenY_centroids,redClusterLabels,greenClusterLabels,colocalisationDistance):
		colocalisedGreenLabels = []
		colocalisedRedLabels = []
		for redX, redY in zip(redX_centroids, redY_centroids):
			distances = np.sqrt(np.square(greenX_centroids-redX) + np.square(greenY_centroids-redY))
			greenMask = (distances<colocalisationDistance)
			if np.sum(greenMask)>0:
				colocalisedGreenLabels.append(np.asarray(greenMask).nonzero())
		for greenX, greenY in zip(greenX_centroids, greenY_centroids):
			distances = np.sqrt(np.square(redX_centroids-greenX) + np.square(redY_centroids-greenY))
			redMask = (distances<colocalisationDistance)
			if np.sum(redMask)>0:
				colocalisedRedLabels.append(np.asarray(redMask).nonzero())
		return colocalisedRedLabels,colocalisedGreenLabels
		
	def plotClusters(dataClustered, clusterLabels, yFOVSize):
		uniqueLabels = set(clusterLabels)
		# Different colour for clusters
		colors = [plt.cm.Spectral(each)
				  for each in np.linspace(0, 1, len(uniqueLabels))]
		random.shuffle(colors)
		for k, col in zip(uniqueLabels, colors):
			if k != -1:
				# Noise not plotted.
				class_member_mask = (clusterLabels == k)
				xy = dataClustered[class_member_mask]	
				
				#yFOVSize used so the image isn't reversed compared to the original image
				plt.plot(xy.iloc[:, 0], yFOVSize - xy.iloc[:, 1], 'o', markerfacecolor=tuple(col), markeredgecolor='k', markersize=5)
		
		n_clusters_ = len(set(clusterLabels)) - (1 if -1 in clusterLabels else 0)
		
		plt.title('Estimated number of clusters: %d' % n_clusters_)
		plt.show()
		return
    
	def hDBScan(dataToCluster, minClusterSize):
	# Entering a minimum cluster size, this will cluster data by the HDBScan algorithm
		
		#Perform the HDBScan
		clusterer = hdbscan.HDBSCAN(minClusterSize, gen_min_span_tree = False)
		clusterer.fit(dataToCluster)
		labels = clusterer.labels_
		
		return labels
		
	def dBScan(dataToCluster, minSamples, epsilon):
	# Epsilon is the distance for something to be considered in the neighbourhood and minSamples is the number of points within epsilon for the DBScan to make a core point
		
		db = DBSCAN(epsilon, minSamples).fit(dataToCluster)
		labels = db.labels_
		
		return labels


"""

				
#Add the results to a dataframe and save them to csv
results = pd.DataFrame({'X position of centroid': x_centroids})
results.insert(1,'Y position of centroid', y_centroids)
results.insert(2,'Area of cluster', areas)
results.insert(3,'Largest diagonal distance', distances)
results.insert(4,'Number of localisations in the cluster', count)
results.to_csv(filepath + 'X31_localisations_centroids_and_areas.csv')




plt.show()
"""
