#extract nodes from Salish Sea Model and average over a year
import math
import os
import csv
import netCDF4 as nc

# function to extract list from csv file
def get_node_list(fileName):
    nodeList = []
    with open(fileName,'r') as read_obj:
        f = csv.reader(read_obj)
        for row in f:
            print(row[0])
            nodeList.append(int(row[0]))
    return nodeList

def main():

    # define file paths
    folder = r"\scratch\05410\eabbott\kelp\SSM"
    fileList = os.listdir(folder)
    nodeList = get_node_list("PugetNodeList.csv") #only nodes within Puget Sound
    resultFileName = r"avgMinDailySal.csv"

    # define variables
    numNodes = len(nodeList)
    numFiles = 364
    badFile = 287 # day 287 file is missing variables
    numHours = 24
    siglay = 0 # surface salinity only? i guess
    FirstSummerDay = 121 # summer days so I can grab Jan - Apr 30 and Nov 1 - Dec 31
    LastSummerDay = 272

    # create a list which contains an empty list for every node
    SalDailyMinList = []
    for i in range(numNodes):
        SalDailyMinList.append([])

    # iterate through every file (day) and calculate
    # the minimum daily salinity at each node
    for iFile in range (numFiles):
        if iFile <= FirstSummerDay or iFile >= LastSummerDay:

            # open the nc file
            print(iFile)
            fn=folder+"\\"+fileList[iFile]
            print(fn)
            ds = nc.Dataset(fn)

            # iterate through nodes of interest
            for iNode in range(numNodes):
                nodeID = nodeList[iNode]

                # find the maximum temperature and put it in the matrix
                SalDailyMin = 0
                for iTime in range(numHours):
                    sal = ds['salinity'][iTime,siglay, nodeID]
                    if sal>SalDailyMin:
                        SalDailyMin = sal
                SalDailyMinList[iNode].append(SalDailyMin)

    # write the results to a csv file
    fResult = open(resultFileName,'w')
    fResult.write("nodeID,avgDailyMinSal\n")
    for iNode in range(numNodes):
        nodeID = nodeList[iNode]
        meandailyminsal = sum(SalDailyMinList[iNode])/len(SalDailyMinList[iNode])
        fResult.write(str(nodeID)+","+str(meandailyminsal)+"\n")

main()
