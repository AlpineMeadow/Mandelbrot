#! /usr/bin/env python3

#A program to create the Mandelbrot set.



def getArgs(parser) :
 
  #Get the parameters
  parser.add_argument('-sX', '--StartX', default = -0.205,
                      help = 'Choose the Starting X coordinate', type = float) 

  parser.add_argument('-sY', '--StartY', default = 1.065,
                      help = 'Choose the Starting Y coordinate', type = float)
  
  parser.add_argument('-eX', '--EndX', default = -0.175,
                      help = 'Choose the Ending X coordinate', type = float)
  
  parser.add_argument('-eY', '--EndY', default = 1.085,
                      help = 'Choose the Ending Y coordinate', type = float) 
  
  parser.add_argument('-nI', '--numIterations', default = 50,
                      help = 'Choose how many iterations are made before stopping.', type = int)

  numGridPointsStr1 = ('Choose the number of grid points. ')
  numGridPointsStr2 = ('This value will give numGridPoints Squared points.')
  numGridPointsStr = numGridPointsStr1 + numGridPointsStr2
  
  parser.add_argument('-nG', '--numGridPoints', default = 100, help = numGridPointsStr,
                      type = int)

  args = parser.parse_args()

  #Generate variables from the inputs.
  numIterations = args.numIterations
  numGridPoints = args.numGridPoints
  startX = args.StartX
  startY = args.StartY
  endX = args.EndX
  endY = args.EndY
  xyLocations = [startX, startY, endX, endY]
  
  return numGridPoints, numIterations, xyLocations
#End of the function getArgs(parser).py

#################################################################################

#################################################################################

def plotMandelbrot(MSet, numIterations, xyLocations, outfile) :
  import numpy as np
  import matplotlib.pyplot as plt
  from matplotlib.backends.backend_pdf import PdfPages
  import matplotlib as mpl
  from math import sqrt
  
  #Set out the coordinates and the set membership indicator.
  x = MSet[:, 0]
  y = MSet[:, 1]
  inSet = MSet[:, 2]
  notInSet = MSet[:, 3]
  
  #generate a color map that only plots the Mandelbrot set itself.  The levels values must not
  #include the values found in inSet.  inSet contains only zeros and ones so the last level
  #must be larger than 1.0 or else this doesn't work.
  inSetLevels = [0.0, 0.25, 1.5]
  inSetColors = ['gold', 'black']
  inSetcmap, inSetnorm = mpl.colors.from_levels_and_colors(inSetLevels, inSetColors)

  #Generate a color map that plots the degree to which the points leave the Mandelbrot set.
  notInSetLevels = list(np.arange(1.01, numIterations + 0.01, 10))
  notInSetColors = ['midnightblue', 'royalblue', 'cornflowerblue', 'lightsteelblue']
  notInSetcmap, notInSetnorm = mpl.colors.from_levels_and_colors(notInSetLevels, notInSetColors)

  #Generate a title string.
  xRes = str(int(sqrt(len(x))))
  yRes = str(int(sqrt(len(y))))  
  titleStr = 'Mandelbrot Set - Resolution : ' + xRes + 'x' + yRes

  #Plot the set.
  xMin = xyLocations[0]
  xMax = xyLocations[2]
  yMin = xyLocations[1]
  yMax = xyLocations[3]
  deltaX = np.abs(xMin - xMax)
  deltaY = np.abs(yMin - yMax)

  #We keep the aspect ratio of the plot the same as the domain.
  figsizeX = 20
  figsizeY = 20*deltaY/deltaX

  plt.figure(figsize = (figsizeX, figsizeY))
  plt.rcParams.update({'font.size': 24})
  ax = plt.axes()

  ax.set_xlim([xMin, xMax])
  ax.set_ylim([yMin, yMax])
  ax.scatter(x, y, marker = '.', c = inSet, cmap = inSetcmap, norm = inSetnorm)
  ax.scatter(x, y, marker = '.', c = notInSet, cmap = notInSetcmap, norm = notInSetnorm)
  plt.plot([-2.0, 2.0], [0.0, 0.0], color = 'black')
  plt.plot([0.0, 0.0], [-2.0, 2.0], color = 'black')
  plt.title(titleStr)

  plt.savefig(outfile)
  #Save the plot to a file.
#  pp = PdfPages(outfile)
#  pp.savefig()
#  pp.close()

  plt.close()
  plt.cla()
  plt.clf()
#End of the function plotMandelbrot.py

######################################################################################

######################################################################################

def getMandelbrotSet(startPoint, numGridPoints, numIterations, xyLocations) :
  import numpy as np
  
  #Allocate an array that will hold the Mandelbrot set.  First column will hold the x component
  #of the point, second column will hold the y component of the point the third column will
  #hold the value zero(not in set) or one(in set) and the fourth column will hold the number of
  #iterations before leaving the set.
  MSet = np.zeros((numGridPoints*numGridPoints, 4))

  #Create a set of coordinates
  x = np.linspace(xyLocations[0], xyLocations[2], numGridPoints)
  y = np.linspace(xyLocations[1], xyLocations[3], numGridPoints)

  #Loop through the gridpoints.  These are the c's in the original formula.
  for i in range(numGridPoints) :
    for j in range(numGridPoints) :

      #Generate the constant and the initial function value.
      c = complex(x[i], y[j])
      z = startPoint**2 + c      
  
      #Iterate the function.  If magnitude of function at location stays bounded, then declare
      #point to be in set, else declare it outside.
      for u in range(0, numIterations) :
        z = z**2 + c
        zMag = abs(z)

        #Check to see if leaving the set.
        if(zMag >= 2.0) :
          inOrOut = 0
          iterationsBeforeLeaving = u
          break
        else :
          inOrOut = 1
          iterationsBeforeLeaving = numIterations
      #End of for loop - for u in range(0, numIterations) :
      
      MSet[numGridPoints*i + j] = x[i], y[j], inOrOut, iterationsBeforeLeaving
    #End of for loop - for j in range(numGridPoints) :
  #End of for loop - for i in range(numGridPoints) ;

  return MSet
#End of the function getMandelbrotSet.py

#####################################################################################

#####################################################################################

#Gather our code in a main() function.
def main() :
  import numpy as np
  import gc
  
  #Set up the original starting point.
  z0 = 0 + 0j

  numGridPoints = 1000
  numIterations = 50

  #Generate x and y coordinates with domain of -2.0 to 2.0 divided by a given number.
  startX = -2.0
  endX = 2.0
  startY = -2.0
  endY = 2.0
  
  numCoordinates = 50
  XCoords = np.linspace(startX, endX, numCoordinates)
  YCoords = np.linspace(startY, endY, numCoordinates)
  xLength = np.linalg.norm(endX - startX)/numCoordinates
  yLength = np.linalg.norm(endY - startY)/numCoordinates

  #Loop through the coordinates.
  for i in range(len(XCoords)) :
    for j in range(len(YCoords)) :

      #Set the x,y locations.
      xyLocations = [XCoords[i], YCoords[j], XCoords[i] + xLength, YCoords[j] + yLength]
      
      #Get the Mandelbrot set.
      MSet = getMandelbrotSet(z0, numGridPoints, numIterations, xyLocations)
  
      #Create a output file name to where the plot will be saved.
      outfilepath = '/home/jdw/Computer/Mandelbrot/Plots/'

      startXStr = ('{:-.3f}'.format(XCoords[i]))
      startYStr = ('{:-.3f}'.format(YCoords[j]))
      endXStr = ('{:-.3f}'.format(XCoords[i] + xLength))
      endYStr = ('{:-.3f}'.format(YCoords[j] + yLength))
      
      filename = ('Mandelbrot_' + startXStr + '_' + endXStr + '_' + startYStr + '_' +
                  endYStr + '_' + str(numIterations) + '.jpg')
      outfile = outfilepath + filename     

      #Now plot the results.
      plotMandelbrot(MSet, numIterations, xyLocations, outfile)

      #I think I am having memory issues.  I will delete the MSet variable and then garbage
      #collect to free up memory.
      del MSet
      gc.collect()

    #End of for loop - for j in range(len(YCoords)) :
  #End of for loop - for i in range(len(XCoords)) :
  
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
  main()
