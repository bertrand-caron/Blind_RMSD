import math
import numpy

def alignPointsOnPoints(point_list1, point_list2):
    cog1 = center_of_mass(point_list1)
    cog2 = center_of_mass(point_list2)
    translate_cog=substract(cog2, cog1)
    print translate_cog
    print "Aligning by translating from center of mass difference"
    return map(lambda point:add(point, translate_cog), point_list1)

def rmsd(point_list1, point_list2):
    rmsd = 0.
    n_points = len(point_list1)
    for point1 in point_list1:
        #print "Point1: {0}".format(point1)
        closest_distance = min([ distance(point1, point2) for point2 in point_list2])
        #print "Closest_distance: {0}".format(closest_distance)
        rmsd += closest_distance
    rmsd = math.sqrt( rmsd / n_points )
    print "RMSD: {0}".format(rmsd)
    return rmsd

def distance(point1, point2):
    return math.sqrt( (point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2 )

def center_of_mass(point_list):
    x = sum([ point[0] for point in point_list])
    y = sum([ point[1] for point in point_list])
    z = sum([ point[2] for point in point_list])
    return map(lambda x:x/len(point_list), [x, y, z])

def substract(point1, point2):
    return [ point1[0]-point2[0], point1[1]-point2[1], point1[2]-point2[2]  ]

def add(point1, point2):
    return [ point1[0]+point2[0], point1[1]+point2[1], point1[2]+point2[2]  ]

if __name__ == "__main__":
    main(point_list1, point_list2)
