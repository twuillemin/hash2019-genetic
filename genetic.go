package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/bits"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"
)

func main() {

	// Init the randomness
	rand.Seed(time.Now().UnixNano())

	processFile("data/a_example")           // Example
	processFile("data/b_lovely_landscapes") // All Horizontal
	processFile("data/c_memorable_moments") // Mixed small
	processFile("data/d_pet_pictures")      // Mixed large
	processFile("data/e_shiny_selfies")     // All vertical

}

// maximumSteps is the maximum number of steps to be run
const maximumSteps = 1000

// populationSize is the size of the population.
const populationSize = 1000

// elitismFactor is the percentage of best parents that are kept in the next generation.
const elitismFactor = 0.1

// useRoulette defines if we use the roulette algorithm for breeding selection. Otherwise use tournament.
const useRoulette = false

// mutationFactor is the factor applied to each gene to determine if it mutates, keep it small
const mutationFactor = 0.0005

// orientation is the orientation of a photo defined as vertical or horizontal.
type orientation int

const (
	horizontal orientation = 1
	vertical   orientation = 2
)

// photoDefinition is the definition of a photo.
type photoDefinition struct {
	index       int
	orientation orientation
	tags        map[string]bool
	tagBits     []uint64
	nbTags      int
}

// scoreResult glues together an index in the population and its score.
type scoreResult struct {
	genomeIndex int
	score       int
}

// childrenResult glues together an index and a children genome
type childrenResult struct {
	childrenId int
	genome     *[]int
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                               MAIN FUNCTION
//
// -------------------------------------------------------------------------------------------------------------------

// processFile processes entirely a file from reading its data to generating the results. This is the heart of the process
func processFile(fileName string) {

	fmt.Printf("Start processing: %v\n", fileName)

	photos, tags := readFile(fileName)

	// Detect if tags are sparse
	sparseTags := false
	if len(tags) > 10*len(*photos) {
		fmt.Printf("Sparse tags detected\n")
		sparseTags = true
	}

	population := createInitialPopulation(len(*photos), populationSize)
	scores := make([]int, populationSize)

	start := time.Now()

	currentMax := -1
	stableMaxCount := 0

	// Obviously there are better strategies than
	for round := 0; round < maximumSteps && stableMaxCount < 10; round++ {

		fmt.Printf("Generation: %v\n", round)

		scores = scorePopulation(photos, population, sparseTags)

		max := analyzeScores(scores)

		// If no better increment the stability counter otherwise reset it
		if max == currentMax {
			stableMaxCount++
		} else {
			currentMax = max
			stableMaxCount = 0
		}

		// Only breed and mutate if stability is not reached
		if stableMaxCount < 10 {

			population = breed(population, scores)

			mutate(population)
		}
	}

	// Rank the scores
	bestScores := make([]scoreResult, len(scores))
	for i := range bestScores {
		bestScores[i] = scoreResult{
			i,
			scores[i],
		}
	}

	sort.Slice(bestScores, func(i, j int) bool {
		return bestScores[i].score > bestScores[j].score
	})

	fmt.Printf("Solution found in %v\n", time.Since(start))

	// Write the results
	writeFile(photos, population[bestScores[0].genomeIndex], fileName+".out")
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                    INITIALIZATION FUNCTIONS
//
// -------------------------------------------------------------------------------------------------------------------

// createInitialPopulation creates an initial population of genomes.
func createInitialPopulation(genomeSize int, populationSize int) []*[]int {

	// Create an array from 0 to len
	base := make([]int, genomeSize)
	for i := 0; i < genomeSize; i++ {
		base[i] = i
	}

	population := make([]*[]int, populationSize)

	// For each genome to generate, shuffle a copy of the base array
	for i := 0; i < populationSize; i++ {
		genome := make([]int, genomeSize)
		perm := rand.Perm(genomeSize)
		for j, v := range perm {
			genome[j] = base[v]
		}
		population[i] = &genome
	}

	return population
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                     FITNESS FUNCTIONS
//
// -------------------------------------------------------------------------------------------------------------------

func scorePopulation(photos *[]photoDefinition, population []*[]int, sparseTags bool) []int {

	scores := make([]int, len(population))

	messages := make(chan scoreResult)

	// Score the population
	start := time.Now()
	for i := 0; i < len(population); i++ {
		go func(genomeIndex int) {
			score := 0
			if sparseTags {
				score = scoreGenomeSparse(photos, population[genomeIndex])
			} else {
				score = scoreGenome(photos, population[genomeIndex])
			}
			messages <- scoreResult{genomeIndex, score}
		}(i)
	}

	// Receive the scores
	for i := 0; i < len(population); i++ {
		result := <-messages
		scores[result.genomeIndex] = result.score
	}

	close(messages)

	fmt.Printf("genome scored in %v\n", time.Since(start))

	return scores
}

// scoreGenome scores a genome. The resulting score "should" be the Google Score. The slideware is generated as a slide
// for each horizontal photo and a slice for each pair of vertical photo. For pair the slide is generated when the
// second member of the pair is found
func scoreGenome(photos *[]photoDefinition, genome *[]int) int {

	// Information about the previous slide (whatever if 2 vertical or a single horizontal)
	previousSlideTags := make([]uint64, len((*photos)[0].tagBits))
	previousSlideNbTags := 0

	// The previous vertical, so that when finding a vertical it is possible to create slide
	var previousVertical *photoDefinition

	// temporary variable for computing a 2 * vertical slide to avoid loop allocation
	verticalSlideTags := make([]uint64, len((*photos)[0].tagBits))
	verticalSlideNbTags := 0

	// The score
	score := 0

	for _, photoIndex := range *genome {

		photo := (*photos)[photoIndex]

		// If initialization of scoring
		if previousSlideTags == nil {
			if photo.orientation == horizontal {
				previousSlideTags = photo.tagBits
			} else {
				if previousVertical == nil {
					previousVertical = &photo
				} else {
					// Sum the current vertical and the previous one to make the tag bits of the fist slide
					previousSlideNbTags = 0
					for i := 0; i < len(previousSlideTags); i++ {
						previousSlideTags[i] = previousVertical.tagBits[i] | photo.tagBits[i]
						previousSlideNbTags += bits.OnesCount64(previousSlideTags[i])
					}
					previousVertical = nil
				}
			}
		} else {
			if photo.orientation == horizontal {
				// Compute the score
				score += getInterestByBits(previousSlideTags, previousSlideNbTags, photo.tagBits, photo.nbTags)
				// Keep the current slide as the previous one
				previousSlideTags = photo.tagBits
				previousSlideNbTags = photo.nbTags
			} else {
				if previousVertical == nil {
					previousVertical = &photo
				} else {
					// Sum the current vertical and the previous one to make the tag bits of the slide
					verticalSlideNbTags = 0
					for i := 0; i < len(previousSlideTags); i++ {
						verticalSlideTags[i] = previousVertical.tagBits[i] | photo.tagBits[i]
						verticalSlideNbTags += bits.OnesCount64(verticalSlideTags[i])
					}
					// Remove the previous vertical
					previousVertical = nil
					// Compute the score
					score += getInterestByBits(previousSlideTags, previousSlideNbTags, verticalSlideTags, verticalSlideNbTags)
					// Keep the current slide as the previous one
					previousSlideTags = verticalSlideTags
					previousSlideNbTags = verticalSlideNbTags
				}
			}
		}
	}

	return score
}

// getInterestByBits computes the score of two slides by their tag vectors
func getInterestByBits(bits1 []uint64, nbBitsOn1 int, bits2 []uint64, nbBitsOn2 int) int {

	commonTags := 0
	for i := 0; i < len(bits1); i++ {
		tmp := bits1[i] & bits2[i]
		commonTags += bits.OnesCount64(tmp)
	}

	minimum := commonTags

	if nbBitSpecific1 := nbBitsOn1 - commonTags; nbBitSpecific1 < minimum {
		minimum = nbBitSpecific1
	}

	if nbBitSpecific2 := nbBitsOn2 - commonTags; nbBitSpecific2 < minimum {
		minimum = nbBitSpecific2
	}

	return minimum
}

// scoreGenome scores a genome. The resulting score "should" be the Google Score. The slideware is generated as a slide
// for each horizontal photo and a slice for each pair of vertical photo. For pair the slide is generated when the
// second member of the pair is found
func scoreGenomeSparse(photos *[]photoDefinition, genome *[]int) int {

	// Information about the previous slide (whatever if 2 vertical or a single horizontal)
	previousSlideTags := make(map[string]bool)

	// The previous vertical, so that when finding a vertical it is possible to create slide
	var previousVertical *photoDefinition

	// The score
	score := 0

	for _, photoIndex := range *genome {

		photo := (*photos)[photoIndex]

		// If initialization of scoring
		if len(previousSlideTags) == 0 {
			if photo.orientation == horizontal {
				previousSlideTags = photo.tags
			} else {
				if previousVertical == nil {
					previousVertical = &photo
				} else {
					// Sum the current vertical and the previous one to make the tag bits of the fist slide
					previousSlideTags = make(map[string]bool)
					for k := range previousVertical.tags {
						previousSlideTags[k] = true
					}
					for k := range photo.tags {
						previousSlideTags[k] = true
					}

					previousVertical = nil
				}
			}
		} else {
			if photo.orientation == horizontal {
				// Compute the score
				score += getInterestByTag(previousSlideTags, photo.tags)

				// Keep the current slide as the previous one
				previousSlideTags = photo.tags
			} else {
				if previousVertical == nil {
					previousVertical = &photo
				} else {
					// Sum the current vertical and the previous one to make the tag bits of the fist slide
					verticalSlideTags := make(map[string]bool)
					for k := range previousVertical.tags {
						verticalSlideTags[k] = true
					}
					for k := range photo.tags {
						verticalSlideTags[k] = true
					}

					// Remove the previous vertical
					previousVertical = nil

					// Compute the score
					score += getInterestByTag(previousSlideTags, verticalSlideTags)

					// Keep the current slide as the previous one
					previousSlideTags = verticalSlideTags
				}
			}
		}
	}

	return score
}

// getInterestByTag computes the score of two slides by their tags
func getInterestByTag(tags1 map[string]bool, tags2 map[string]bool) int {

	commonTags := 0
	for tag := range tags1 {
		if _, ok := tags2[tag]; ok {
			commonTags++
		}
	}

	minimum := commonTags

	if nbTagSpecific1 := len(tags1) - commonTags; nbTagSpecific1 < minimum {
		minimum = nbTagSpecific1
	}

	if nbTagSpecific2 := len(tags2) - commonTags; nbTagSpecific2 < minimum {
		minimum = nbTagSpecific2
	}

	return minimum
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                    BREEDING AND CROSS OVER FUNCTIONS
//
// -------------------------------------------------------------------------------------------------------------------

// breed takes a population and its score and generate a new population to be scored. The breed is elitist and keeps
// the top `elitismFactor` of the best parents.
func breed(population []*[]int, scores []int) []*[]int {

	start := time.Now()

	// Rank the scores
	bestScores := make([]scoreResult, len(scores))
	for i := range bestScores {
		bestScores[i] = scoreResult{
			i,
			scores[i],
		}
	}

	sort.Slice(bestScores, func(i, j int) bool {
		return bestScores[i].score > bestScores[j].score
	})

	results := make([]*[]int, len(population))

	elitismKept := int(float64(len(scores)) * elitismFactor)

	// Keep good parents
	for i := 0; i < elitismKept; i++ {
		results[i] = population[bestScores[i].genomeIndex]
	}

	messages := make(chan childrenResult)

	// Create new children
	for i := elitismKept; i < len(population); i++ {
		go func(childrenId int) {
			children := fornicate(population, scores)
			messages <- childrenResult{childrenId, children}
		}(i)
	}

	// Receive the children
	for i := elitismKept; i < len(population); i++ {
		result := <-messages
		results[result.childrenId] = result.genome
	}

	close(messages)

	fmt.Printf("population breeded in %v\n", time.Since(start))

	return results
}

// fornicate chooses two parents and make a new children
func fornicate(population []*[]int, scores []int) *[]int {

	normalizedScores := normalizeScore(scores)

	var breeder1, breeder2 int
	if useRoulette {
		breeder1 = selectBreederByRoulette(normalizedScores, -1)
		breeder2 = selectBreederByRoulette(normalizedScores, breeder1)
	} else {
		breeder1 = selectBreederByTournament(normalizedScores, -1)
		breeder2 = selectBreederByTournament(normalizedScores, breeder1)
	}
	return crossOver(population[breeder1], population[breeder2])
}

func normalizeScore(scores []int) []float64 {

	// Find the max and min of the scores
	max := -1
	min := 10000000
	for i := range scores {
		if scores[i] > max {
			max = scores[i]
		}
		if scores[i] < min {
			min = scores[i]
		}
	}

	rangeScore := max - min

	// Normalize the score
	normalized := make([]float64, len(scores))
	sumNormalized := 0.0
	for i := range scores {
		// normalizedScore is from 0 to 1
		normalizedScore := float64(scores[i]-min) / float64(rangeScore)
		normalized[i] = normalizedScore * normalizedScore
		sumNormalized += normalized[i]
	}

	return normalized
}

// selectBreederByRoulette selects randomly a breeder by using a roulette algorithm.
func selectBreederByRoulette(scores []float64, forbiddenBreeder int) int {

	sumScores := 0.0
	for _, v := range scores {
		sumScores += v
	}

	// Rand limits are inclusive
	random := rand.Float64() * sumScores
	sum := 0.0

	// While the selected breeder is the forbidden bridder, find a breeder
	breederSelected := forbiddenBreeder
	for breederSelected == forbiddenBreeder {
		for i, v := range scores {
			sum += v
			if sum >= random {
				// A breeder was found, quit the search loop
				breederSelected = i
				break
			}
		}
	}

	// Return the selected breeder
	return breederSelected
}

// selectBreederByRoulette selects randomly a breeder by using a tournament algorithm.
func selectBreederByTournament(scores []float64, forbiddenBreeder int) int {

	tournamentSize := len(scores) / 8

	bestScore := -1.0
	var bestBreeder int

	for i := 0; i < tournamentSize; i++ {
		possibleBreeder := rand.Intn(len(scores))
		if scores[possibleBreeder] > bestScore && possibleBreeder != forbiddenBreeder {
			bestScore = scores[possibleBreeder]
			bestBreeder = possibleBreeder
		}
	}

	return bestBreeder
}

// crossOver crosses two parents for generating a new children.
func crossOver(parent1 *[]int, parent2 *[]int) *[]int {

	size := len(*parent1)

	start := rand.Intn(size)
	length := rand.Intn(size)
	if start+length > size {
		length = size - start
	}

	children := make([]int, size)
	alreadyFromParent1 := make(map[int]bool, length)

	for i := start; i < start+length; i++ {
		parent1Value := (*parent1)[i]
		children[i] = parent1Value
		alreadyFromParent1[parent1Value] = true
	}

	indexWritingParent2 := 0
	for i := 0; i < size; i++ {
		parent2Value := (*parent2)[i]

		// If already handled by parent1, do not duplicate
		if _, ok := alreadyFromParent1[parent2Value]; ok {
			continue
		}

		// If want to write in the zone filled by parent1, jump the writing at the end
		if start <= indexWritingParent2 && indexWritingParent2 < start+length {
			indexWritingParent2 = start + length
		}

		// If the writing is to be done after the end of the children, this implies that all
		// items remaining in parent2 were already written by parent1
		if indexWritingParent2 >= size {
			break
		}

		children[indexWritingParent2] = parent2Value
		indexWritingParent2++
	}

	return &children
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                    MUTATION FUNCTIONS
//
// -------------------------------------------------------------------------------------------------------------------

// mutate applies random mutation on a population. Mutations are made by swapping the mutated gene with another.
func mutate(population []*[]int) {

	start := time.Now()

	// Parallel mutations are way slower (3*times) than brute force mutation
	// So just brute force it
	for _, genome := range population {
		mutateGenome(genome)
	}

	fmt.Printf("population mutated in %v\n", time.Since(start))
}

// mutate applies random mutation on a population. Mutations are made by swapping the mutated gene with another.
func mutateGenome(genome *[]int) {

	genomeSize := len(*genome)
	// For all genes
	for i := range *genome {
		// If the gene is to be mutated
		if rand.Float64() < mutationFactor {

			// Make a random swap, but not with self
			if idxTarget := rand.Intn(genomeSize); idxTarget != i {
				temp := (*genome)[idxTarget]
				(*genome)[idxTarget] = (*genome)[i]
				(*genome)[i] = temp
			}
		}
	}
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                    INPUT AND OUTPUT FUNCTIONS
//
// -------------------------------------------------------------------------------------------------------------------

// readFiles reads the definition of the photos
func readFile(fileName string) (*[]photoDefinition, map[string]int) {

	// Get the tags
	allTags := readTagsFromFile(fileName)
	// Number of uint64 needed by image
	vectorTagSize := int(math.Ceil(float64(len(allTags)) / 64.0))

	fmt.Printf("Found %v distinct tags (vector size is %v)\n", len(allTags), vectorTagSize)

	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}

	defer func() {
		if e := file.Close(); e != nil {
			log.Panic(err)
		}
	}()

	scanner := bufio.NewScanner(file)

	// Read the header line
	scanner.Scan()
	headerLine := scanner.Text()

	nbImages, err := strconv.Atoi(headerLine)
	if err != nil {
		log.Fatal(err)
	}

	photos := make([]photoDefinition, nbImages)

	for i := 0; i < nbImages; i++ {
		scanner.Scan()
		line := scanner.Text()

		if e := scanner.Err(); e != nil {
			log.Fatal(e)
		}

		lineItems := strings.Split(line, " ")

		photoOrientation := horizontal
		if lineItems[0] != "H" {
			photoOrientation = vertical
		}

		nbTags, e := strconv.Atoi(lineItems[1])
		if e != nil {
			log.Fatal(e)
		}

		tags := make(map[string]bool, nbTags)
		vectorTags := make([]uint64, vectorTagSize)

		for j := 0; j < nbTags; j++ {

			tag := lineItems[j+2]

			tags[tag] = true

			tagIndex := allTags[tag]

			// Compute the position of the bit
			uint64idx := uint(tagIndex) >> 6
			posInUint64 := uint(tagIndex % 64)
			bitToSet := uint64(1) << posInUint64

			// Set the bit
			vectorTags[uint64idx] |= bitToSet
		}

		photos[i] = photoDefinition{
			index:       i,
			orientation: photoOrientation,
			tags:        tags,
			tagBits:     vectorTags,
			nbTags:      len(tags),
		}
	}

	fmt.Printf("Read %v images\n", len(photos))

	return &photos, allTags
}

// readFiles reads the definition of the photos
func readTagsFromFile(fileName string) map[string]int {

	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}

	defer func() {
		if e := file.Close(); e != nil {
			log.Panic(err)
		}
	}()

	scanner := bufio.NewScanner(file)

	// Read the header line
	scanner.Scan()
	headerLine := scanner.Text()

	nbImages, err := strconv.Atoi(headerLine)
	if err != nil {
		log.Fatal(err)
	}

	tags := make(map[string]int)

	for i := 0; i < nbImages; i++ {
		scanner.Scan()
		line := scanner.Text()

		if e := scanner.Err(); e != nil {
			log.Fatal(e)
		}

		lineItems := strings.Split(line, " ")

		for tagIdx := 2; tagIdx < len(lineItems); tagIdx++ {
			tag := lineItems[tagIdx]
			if _, ok := tags[tag]; !ok {
				idx := len(tags)
				tags[tag] = idx
			}
		}
	}

	return tags
}

// writeFile writes a genome to an output file
func writeFile(photos *[]photoDefinition, genome *[]int, fileName string) {

	previousVertical := -1

	// Make an array for keeping the lines
	lines := make([]string, 0)
	lines = append(lines, "count")

	for _, photoIndex := range *genome {

		photo := (*photos)[photoIndex]

		if photo.orientation == horizontal {
			lines = append(lines, fmt.Sprintf("%v\n", photoIndex))
		} else {
			if previousVertical == -1 {
				previousVertical = photoIndex
			} else {
				lines = append(lines, fmt.Sprintf("%v %v\n", previousVertical, photoIndex))
				previousVertical = -1
			}
		}
	}

	// Adjust the size
	lines[0] = fmt.Sprintf("%v\n", len(lines)-1)

	file, err := os.Create(fileName)
	if err != nil {
		log.Panic(err)
	}

	defer func() {
		if e := file.Close(); e != nil {
			log.Panic(err)
		}
	}()

	for _, line := range lines {
		_, err = file.WriteString(line)
		if err != nil {
			log.Panic(err)
		}
	}
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                           UTILITIES
//
// -------------------------------------------------------------------------------------------------------------------

// analyzeScores displays information about the score computed (mean, avg, etc.) and return the best score
func analyzeScores(scores []int) int {

	min := 1000000
	max := -1
	sum := 0

	for _, v := range scores {
		if v < min {
			min = v
		}
		if v > max {
			max = v
		}
		sum += v
	}

	mean := float64(sum) / float64(len(scores))

	sumDiffMeanSquared := 0.0
	for _, v := range scores {
		diff := float64(v) - mean
		sumDiffMeanSquared += diff * diff
	}

	deviation := math.Sqrt(sumDiffMeanSquared / float64(len(scores)))

	fmt.Printf("Min: %v, Max: %v, Mean: %v, Deviation: %v\n", min, max, mean, deviation)

	return max
}
