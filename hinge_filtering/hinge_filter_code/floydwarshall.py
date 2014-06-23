#!/usr/bin/env python

#ryan g. coleman ryangc@mail.med.upenn.edu
# floyd warshall code
#adapted from CLRS of course
#O(n^3) all pairs shortest paths
#just gives distances currently, not actual paths (no Pi matrix)

import unittest
import sys  # for maxsize

def makeMatrix(size, infinity=sys.maxsize):
  if size > 0:
    retMat = []
    for count in xrange(size):
      oneRow = [infinity for count in xrange(size)]
      retMat.append(oneRow)
    return retMat
  else:
    return False

def floydWarshall(neighbors, infinity=sys.maxsize):
  '''
  this version takes a dictionary of neighbors and distances. format is:
  startnode->[[neighbor, dist], [neighbor, dist], [...]]
  '''
  size = len(neighbors)
  oldMat = makeMatrix(size, infinity)
  orderKeys = {}
  #now initialize from the neighbors
  orderedKeys = neighbors.keys()
  orderedKeys.sort()
  for order, key in enumerate(orderedKeys):
    orderKeys[key] = order
  for diagonal in xrange(size):
    oldMat[diagonal][diagonal] = 0
  for key in orderedKeys:
    neighborList = neighbors[key]
    for neigh, dist in neighborList:
      oldMat[orderKeys[key]][orderKeys[neigh]] = dist
      #symmetric case be done later
  newMat = oldMat[:]
  for mac in xrange(size):  # this is the main loop of the code
    oldMat = newMat[:]      # copy the old matrix, work with the new one
    for row in xrange(size):
      for col in xrange(size):
        newMat[row][col] = min(   # update all connections 1 away from previous
            oldMat[row][col], oldMat[row][mac] + oldMat[mac][col])
  return newMat, orderKeys

class TestFloydWarshall(unittest.TestCase):
  '''unit tests of floyd warshall code'''

  def setUp(self):
    self.neighborTest1 = {}
    self.neighborTest1[0] = [[1, 2]]
    self.neighborTest1[1] = [[0, 2], [2, 1]]
    self.neighborTest1[2] = [[1, 1], [3, 4]]
    self.neighborTest1[3] = [[4, 5], [2, 4]]
    self.neighborTest1[4] = [[3, 5]]
    self.neighborTest1[5] = [[6, 2]]
    self.neighborTest1[6] = [[5, 2]]
    #more simple
    self.neighborTest2 = {}
    self.neighborTest2[0] = [[1, 2]]
    self.neighborTest2[1] = [[0, 2], [2, 1]]
    self.neighborTest2[2] = [[1, 1]]
    #empty test
    self.neighborTest3 = {}
    self.neighborTest3[0] = []
    self.neighborTest3[1] = []
    #float test
    self.neighborTest4 = {}
    self.neighborTest4[0] = [[1, 2.5]]
    self.neighborTest4[1] = [[0, 2.5], [2, 1.2]]
    self.neighborTest4[2] = [[1, 1.2]]

  def test_1(self):
    matrix, orderKeys = floydWarshall(self.neighborTest1)
    self.assertEqual(matrix, [
        [0, 2, 3, 7, 12, sys.maxsize, sys.maxsize],
        [2, 0, 1, 5, 10, sys.maxsize, sys.maxsize],
        [3, 1, 0, 4, 9, sys.maxsize, sys.maxsize],
        [7, 5, 4, 0, 5, sys.maxsize, sys.maxsize],
        [12, 10, 9, 5, 0, sys.maxsize, sys.maxsize],
        [
            sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize,
            0, 2],
        [
            sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize,
            2, 0]])
    self.assertEqual(
        sorted(list(orderKeys.iteritems())),
        sorted(list({0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}.iteritems())))

  def test_2(self):
    matrix, orderKeys = floydWarshall(self.neighborTest2)
    self.assertTrue(matrix == [
        [0, 2, 3], [2, 0, 1], [3, 1, 0]])
    self.assertEqual(
        sorted(list(orderKeys.iteritems())),
        sorted(list({0: 0, 1: 1, 2: 2}.iteritems())))

  def test_3(self):
    matrix, orderKeys = floydWarshall(self.neighborTest3)
    self.assertTrue(matrix == [
        [0, sys.maxsize], [sys.maxsize, 0]])
    self.assertEqual(
        sorted(list(orderKeys.iteritems())),
        sorted(list({0: 0, 1: 1}.iteritems())))

  def test_4(self):
    matrix, orderKeys = floydWarshall(self.neighborTest4)
    self.assertTrue(matrix == [
        [0, 2.5, 3.7], [2.5, 0, 1.2], [3.7, 1.2, 0]])
    self.assertEqual(
        sorted(list(orderKeys.iteritems())),
        sorted(list({0: 0, 1: 1, 2: 2}.iteritems())))

if __name__ == '__main__':  # if run, just do tests
  suite = unittest.TestLoader().loadTestsFromTestCase(TestFloydWarshall)
  unittest.TextTestRunner(verbosity=9).run(suite)
