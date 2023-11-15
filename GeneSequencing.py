import time
import math
import random

MATCH = -3
SUB = 1
INDEL = 5

TOP = "TOP"
DIAGONAL = "DIAGONAL"
LEFT = "LEFT"

MAXINDELS = 3

class GeneSequencing: 
	def align(self, sequence1, sequence2, banded, alignmentLength):
		# Set banded and maximum characters for alignment.
		self.banded = banded
		self.MaxCharactersToAlign = alignmentLength

		# Choose alignment method based on 'banded' flag.
		if banded: 
			# Perform banded alignment and retrieve results.
			alignmentCost, alignment1, alignment2 = self.banded_alignment(sequence1, sequence2, alignmentLength)
		else:
			# Perform unrestricted alignment and retrieve results.
			alignmentCost, alignment1, alignment2 = self.unrestrictedAlignment(sequence1, sequence2, alignmentLength)
		
		# Return alignment information with a limited length of 100 characters.
		return {'align_cost': alignmentCost, 'seqi_first100': '{}'.format(alignment1[:100]), 'seqj_first100': '{}'.format(alignment2[:100])}
	
	def unrestrictedAlignment(self, sequence1, sequence2, alignmentLength):
		# Initialize sequences with space and matrices
		sequence1 = " " + sequence1[:alignmentLength]
		sequence2 = " " + sequence2[:alignmentLength]
		
		# Initialize scoring and path matrices
		scoringMatrix = [[None for _ in range(len(sequence2))] for _ in range(len(sequence1))]
		pathMatrix = [[None for _ in range(len(sequence2))] for _ in range(len(sequence1))]

		# Set base case values for scoring and path matrices
		scoringMatrix[0][0] = 0
		pathMatrix[0][0] = None

		# Fill the first row of scoring and path matrices
		for col in range(1, len(sequence2)):
			scoringMatrix[0][col] = INDEL * col
			pathMatrix[0][col] = LEFT  # Left path for indel

		# Fill the first column of scoring and path matrices
		for row in range(1, len(sequence1)):
			scoringMatrix[row][0] = INDEL * row
			pathMatrix[row][0] = TOP  # Top path for indel

		# Populate scoring and path matrices
		for row in range(1, len(sequence1)):
			for col in range(1, len(sequence2)):
				# Time Complexity for this part: O(1)
				alignmentCost = scoringMatrix[row][col-1] + INDEL
				pathMatrix[row][col] = LEFT  # Left path for indel

				# Check and update for top path
				# Time Complexity for this part: O(1)
				if scoringMatrix[row-1][col] + INDEL < alignmentCost:
					alignmentCost = scoringMatrix[row-1][col] + INDEL
					pathMatrix[row][col] = TOP  # Top path for indel

				# Check and update for diagonal path
				# Time Complexity for this part: O(1)
				if scoringMatrix[row-1][col-1] + self.match(row, col, sequence1, sequence2) < alignmentCost:
					alignmentCost = scoringMatrix[row-1][col-1] + self.match(row, col, sequence1, sequence2)
					pathMatrix[row][col] = DIAGONAL  # Diagonal path for match

				scoringMatrix[row][col] = alignmentCost

		# Initialize pointers for traceback
		row = len(sequence1) - 1
		col = len(sequence2) - 1

		alignedSequence1 = sequence1
		alignedSequence2 = sequence2

		# Initialize alignment cost
		alignmentCost = scoringMatrix[len(sequence1)-1][len(sequence2)-1]

		# Traceback to construct aligned sequences
		while True:
			# Time Complexity for this part: O(len(sequence1) + len(sequence2))
			if row != 0 and col != 0 and pathMatrix is None:
				break
			elif row == 0 and col == 0:
				break
			else:
				# Diagonal path
				# Time Complexity for this part: O(1)
				if pathMatrix[row][col] == DIAGONAL:
					row -= 1
					col -= 1
					continue
				# Left path
				# Time Complexity for this part: O(1)
				elif pathMatrix[row][col] == LEFT:
					alignedSequence1 = alignedSequence1[:row + 1] + "-" + alignedSequence1[row + 1:]
					col -= 1
					continue
				# Top path
				# Time Complexity for this part: O(1)
				elif pathMatrix[row][col] == TOP:
					alignedSequence2 = alignedSequence2[:col + 1] + "-" + alignedSequence2[col + 1:]
					row -= 1

		# Return alignment cost and aligned sequences
		# Time Complexity for this part: O(1)
		return alignmentCost, alignedSequence1[1:], alignedSequence2[1:]

	def match(self, sequence1Index, sequence2Index, sequence1, sequence2):
		# Calculate match or substitution score based on whether characters at given indices are equal | Time Complexity: O(1)	
		return MATCH if sequence1[sequence1Index] == sequence2[sequence2Index] else SUB

	def bandedAlignment(self, sequence1, sequence2, alignmentLength):
		# Determine effective lengths for sequences
		rowLength = min(alignmentLength + 1, len(sequence1) + 1)  # O(1)
		colLength = min(alignmentLength + 1, len(sequence2) + 1)  # O(1)

		# Truncate sequences based on effective lengths
		sequence1 = sequence1[:rowLength - 1]  # O(rowLength)
		sequence2 = sequence2[:colLength - 1]  # O(colLength)

		# Check for significant length difference, return if so
		if abs(rowLength - colLength) > 3:  # O(1)
			return float('inf'), '', ''

		# Set band parameters
		bandWidth = 3  # O(1)
		bandOffset = 7  # O(1)

		# Initialize scoring and path matrices
		scoringMatrix = [[None for _ in range(bandOffset)] for _ in range(rowLength)]  # O(rowLength * bandOffset)
		pathMatrix = [[None for _ in range(bandOffset)] for _ in range(rowLength)]  # O(rowLength * bandOffset)
		scoringMatrix[0][0] = 0  # O(1)

		# Populate scoring and path matrices
		for row in range(rowLength):  # O(rowLength * bandOffset)
			# Define offset for sequence2 based on current row
			sequence2Offset = ' ' + sequence2[:bandOffset] if row < 4 else sequence2[row - bandWidth - 1:bandOffset + (row - bandWidth)]  # O(bandOffset)

			for col in range(bandOffset):  # O(bandOffset)
				# Define offset for sequence1 based on current row
				sequence1Offset = sequence1[row - 1]  # O(1)
				sequence2Col = sequence2Offset[col]  # O(1)

				# Initialize scores for left, top, and diagonal paths
				left, top, diagonal = 0, 0, 0  # O(1)

				# Update scores based on position in matrix
				if row < 4:  # O(1)
					left = scoringMatrix[row][col - 1]  # O(1)
					top = scoringMatrix[row - 1][col]  # O(1)
					diagonal = scoringMatrix[row - 1][col - 1]  # O(1)
				elif col == 6:  # O(1)
					left = scoringMatrix[row][col - 1]  # O(1)
					top = float('inf')  # O(1)
					diagonal = scoringMatrix[row - 1][col]  # O(1)
				elif col == 0:  # O(1)
					left = float('inf')  # O(1)
					top = scoringMatrix[row - 1][col + 1]  # O(1)
					diagonal = scoringMatrix[row - 1][col]  # O(1)
				else:  # O(1)
					left = scoringMatrix[row][col - 1]  # O(1)
					top = scoringMatrix[row - 1][col + 1]  # O(1)
					diagonal = scoringMatrix[row - 1][col]  # O(1)

				# Update diagonal score based on match or substitution
				diagonal = diagonal + SUB if sequence1Offset != sequence2Col else diagonal + MATCH  # O(1)

				# Update scoring matrix and path matrix
				if left + INDEL <= diagonal and left + INDEL <= top + INDEL:  # O(1)
					scoringMatrix[row][col] = left + INDEL  # O(1)
					pathMatrix[row][col] = LEFT  # O(1)
				elif top + INDEL <= diagonal and top + INDEL <= left + INDEL:  # O(1)
					scoringMatrix[row][col] = top + INDEL  # O(1)
					pathMatrix[row][col] = TOP  # O(1)
				else:  # O(1)
					scoringMatrix[row][col] = diagonal  # O(1)
					pathMatrix[row][col] = DIAGONAL  # O(1)

		# Extract the final score from the last row
		for score in scoringMatrix[-1]:  # O(bandOffset)
			if score is not None:  # O(1)
				finalScore = score  # O(1)

		# Reconstruct the alignments using path matrix
		alignment1, alignment2 = self.makeAlignment(pathMatrix, sequence1, sequence2, alignmentLength)  # O(rowLength)

		# Return final score and aligned sequences
		return finalScore, alignment1[1:], alignment2[1:]  # O(1)

	def makeAlignment(self, pathMatrix, sequence1, sequence2, alignmentLength):
		# Set descriptive variable for the threshold
		threshold = 3  # O(1)

		# Initialize row and column indices
		row = len(sequence1) - 1  # O(1)
		col = 0  # O(1)

		# Initialize alignment sequences
		sequence1Alignment = sequence1[:alignmentLength]  # O(alignmentLength)
		sequence2Alignment = sequence2[:alignmentLength]  # O(alignmentLength)

		# Find the end point in the path matrix
		while pathMatrix[row][col] is not None:  # O(alignmentLength)
			col += 1  # O(1)

		col -= 1  # O(1)
		pathMatrix[0][0] = 'START'  # O(1)
		path = pathMatrix[row][col]  # O(1)

		# Reconstruct alignments based on the path
		while path != 'START':  # O(alignmentLength)
			if path == DIAGONAL:  # O(1)
				if row > threshold:  
					row -= 1  
					path = pathMatrix[row][col]  
				else:  # O(1)
					row -= 1  
					col -= 1  
					path = pathMatrix[row][col]  
			elif path == LEFT:  # O(1)
				col -= 1  
				path = pathMatrix[row][col]  
				sequence1Alignment = sequence1Alignment[:row] + "-" + sequence1Alignment[row:]  # 
			elif path == TOP:  # O(1)
				if row > threshold:  
					row -= 1  
					col += 1  
					path = pathMatrix[row][col]  
				else:  
					row -= 1  
					path = pathMatrix[row][col]  

				if row > threshold + 1:  
					sequence2Alignment = sequence2Alignment[:row - threshold + col] + "-" + sequence2Alignment[row - threshold + col:]   
				else:  
					sequence2Alignment = sequence2Alignment[:row - 1] + "-" + sequence2Alignment[row - 1:]  
			else:  # O(1)
				path = "START"  

		return sequence1Alignment, sequence2Alignment  # O(alignmentLength)