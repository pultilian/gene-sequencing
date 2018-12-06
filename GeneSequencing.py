#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF

else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		self.alignA = ""
		self.alignB = ""
		pass

	
# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		results = []
		print(align_length)

		for i in range(len(sequences)):
			jresults = []
			for j in range(len(sequences)):

				if(j < i):
					s = {}
				else:
					score = 0
					grid = None
					lengthY = len(sequences[i]) - 1
					lengthX = len(sequences[j]) - 1
					if lengthX > align_length-2:
						lengthX = align_length
					if lengthY > align_length-2:
						lengthY = align_length
					if banded:
						grid = self.calculateSequenceBound(sequences[i][0:align_length],sequences[j][0:align_length])
						if grid == None:
							score = math.inf
						else:
							score = grid[-1][-1].value
					else:
						grid = self.calculateSequenceUnBound(sequences[i][0:align_length], sequences[j][0:align_length])
						if grid == None:
							score = math.inf
						else:
							score = grid[-1][-1].value
					if grid == None:
						self.alignA = "None"
						self.alignB = "None"
					else:
						self.createAlignment(grid, sequences[j][0:align_length], sequences[i][0:align_length], banded)
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2


					alignment1 = self.alignA[0:100] + ''.format(i+1,
						len(sequences[i]), align_length, ',BANDED' if banded else '')
					print("----------------------------------")
					print(i)
					print(self.alignA[0:100])
					alignment2 = self.alignB[0:100] + ''.format(j+1,
						len(sequences[j]), align_length, ',BANDED' if banded else '')
					print(j)
					print(self.alignB[0:100])
					print("----------------------------------")
###################################################################################################					
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
					table.repaint()	
				jresults.append(s)
			results.append(jresults)
			# self.calculateSequenceBound(sequences[0], sequences[0])
		return results


	def createAlignment(self, grid, StringA, StringB, banded):
		lengthA = len(StringA)
		lengthB = len(StringB)
		breakLoop = False
		alignA = ""
		alignB = ""

		# used for banded calculations
		symbolicX = lengthB

		if banded:
			current = grid[-1][-1]
		else:
			current = grid[lengthA][lengthB]
		while not breakLoop:
			curA = StringA[current.y-1]
			curB = StringB[current.x-1]

			# the end of the path
			if current.parentx == -1 and current.parenty == -1:
				break
			# match and sub calculations
			if current.type == "match" or current.type == "sub":
				alignA += curA
				if banded:
					alignB += StringB[symbolicX - 1]
					symbolicX += -1
				else:
					alignB += curB
			# insert and deletion up
			elif current.type == "indelUp":
				alignA += curA
				alignB += "-"
			# insert and deletion left
			elif current.type == "indelLeft":
				alignA += "-"
				if banded:
					alignB += StringB[symbolicX - 1]
					symbolicX += -1
				else:
					alignB += curB
			current = grid[current.parenty][current.parentx]
		# reverse the string
		self.alignA = alignA[::-1]
		self.alignB = alignB[::-1]

	def calculateSequenceBound(self, StringA, StringB):
		boxList = []
		y = 0
		x = 0
		breakLoop = False
		# loop for the y axis
		while not breakLoop:
			# setup and edge cases
			endIndex = 7 + y - 4
			beginIndex = endIndex - 6
			if(beginIndex < 0):
				beginIndex = 0
			if endIndex > len(StringA):
				endIndex = len(StringA)
			if endIndex == beginIndex:
				break
			elif y > len(StringB):
				break
			else:
				rowList = []
				curIndex = beginIndex
				x = 0

				# Row Operations
				while curIndex <= endIndex:
					canIndelUp = True
					canIndelLeft = True
					canDiag = True
					compareUp = math.inf
					compareDiag = math.inf
					compareLeft = math.inf

					# first
					if x == 0 and y == 0:
						rowList.append(Box(0,-1,-1,0,0,"match"))
					else:
						# Top Row
						if y == 0:
							canIndelUp = False
							canDiag = False

						# in first section
						if beginIndex == 0:
							# Left most Value
							if x == 0:
								canIndelLeft = False
								canDiag = False

							# end of row
							if curIndex == endIndex:
								canIndelUp = False


							# conditions for operations
							if canIndelLeft:
								compareLeft = rowList[x-1].value + 5
							if canIndelUp:
								compareUp = boxList[y - 1][x].value + 5
							diagString = ""
							if canDiag:
								if StringA[curIndex - 1] == StringB[y - 1]:
									compareDiag = boxList[y - 1][x-1].value - 3
									diagString = "match"
								else:
									compareDiag = boxList[y - 1][x-1].value + 1
									diagString = "sub"

							# compareUp
							if compareUp <= compareDiag and compareUp <= compareLeft:
								rowList.append(Box(compareUp, x, y - 1, x, y, "indelUp"))

							# compareLeft
							elif compareLeft <= compareDiag and compareLeft <= compareUp:
								rowList.append(Box(compareLeft, x - 1, y, x, y, "indelLeft"))

							# compareDiag
							elif compareDiag <= compareLeft and compareDiag <= compareUp:
								rowList.append(Box(compareDiag, x - 1, y - 1, x, y, diagString))

						# Not in first set
						else:
							# edge cases
							if curIndex == endIndex:
								canIndelUp = False
							if curIndex == beginIndex:
								canIndelLeft = False

							# checks operations
							if canIndelLeft:
								compareLeft = rowList[x-1].value + 5
							if canIndelUp:
								compareUp = boxList[y - 1][x+1].value + 5
							diagString = ""
							compareDiag = math.inf
							if canDiag:
								if StringA[curIndex - 1] == StringB[y - 1]:
									compareDiag = boxList[y - 1][x].value - 3
									diagString = "match"
								else:
									compareDiag = boxList[y - 1][x].value + 1
									diagString = "sub"

							# compareUp
							if compareUp <= compareDiag and compareUp <= compareLeft:
								rowList.append(Box(compareUp, x+1, y - 1, x, y, "indelUp"))

							# compareLeft
							elif compareLeft <= compareDiag and compareLeft <= compareUp:
								rowList.append(Box(compareLeft, x - 1, y, x, y, "indelLeft"))

							# compareDiag
							elif compareDiag <= compareLeft and compareDiag <= compareUp:
								rowList.append(Box(compareDiag, x, y - 1, x, y, diagString))

					curIndex += 1
					x += 1
				y += 1
				boxList.append(rowList)
		if y < len(StringB):
			boxList = None
		return boxList


	def calculateSequenceUnBound(self, StringA, StringB):
		values = []

		for y in range(len(StringB)+1):
			rowVal = []
			for x in range(len(StringA)+1):
				# end case
				if x == 0 and y == 0:
					rowVal.append(Box(0,-1,-1, 0, 0, "match"))
					continue
				# edge cases
				if x == 0:
					rowVal.append(Box(y*5,x,y-1,x,y, "indelUp"))
				elif y == 0:
					rowVal.append(Box(5*x,x-1,y,x,y, "indelLeft"))
				else:
					# calculate the various values
					diagString = ""
					compareUp = values[y-1][x].value + 5
					compareLeft = rowVal[x-1].value + 5
					compareDiag = math.inf
					if StringA[x-1] == StringB[y-1]:
						compareDiag = values[y-1][x-1].value - 3
						diagString = "match"
					else:
						compareDiag = values[y - 1][x - 1].value + 1
						diagString = "sub"

					# compareUp
					if compareUp <= compareDiag and compareUp <= compareLeft:
						rowVal.append(Box(compareUp,x,y-1,x,y, "indelUp"))

					# compareLeft
					elif compareLeft <= compareDiag and compareLeft <= compareUp:
						rowVal.append(Box(compareLeft,x-1,y,x,y, "indelLeft"))

					# compareDiag
					elif compareDiag <= compareLeft and compareDiag <= compareUp:
						rowVal.append(Box(compareDiag, x - 1, y - 1,x,y, diagString))
			values.append(rowVal)
		return values

class Box:
	def __init__(self, value, parentx, parenty, x, y, type):
		# Value of the box
		self.value = value
		# x and y coordinates of its parent. a negative value means it has no parents
		self.parentx = parentx
		self.parenty = parenty
		self.x = x
		self.y = y
		self.type = type

