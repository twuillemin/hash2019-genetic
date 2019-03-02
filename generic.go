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

	//processFile("data/a_example") // Example
	//processFile("data/b_lovely_landscapes") // All Horizontal
	//processFile("data/c_memorable_moments") // Mixed small
	processFile("data/d_pet_pictures") // Mixed large
	//processFile("data/e_shiny_selfies") // All vertical

}

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

	photos := readFile(fileName)
	fmt.Printf("Read %v images\n", len(*photos))

	population := createInitialPopulation(len(*photos), populationSize)
	scores := make([]int, populationSize)

	for round := 0; round < 100; round++ {

		fmt.Printf("Generation: %v\n", round)

		scores = scorePopulation(photos, population)

		analyzeScores(scores)

		population = breed(population, scores)

		mutate(population)
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

func scorePopulation(photos *[]photoDefinition, population []*[]int) []int {

	scores := make([]int, len(population))

	messages := make(chan scoreResult)

	// Score the population
	start := time.Now()
	for i := 0; i < len(population); i++ {
		go func(genomeIndex int) {
			score := scoreGenome(photos, population[genomeIndex])
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
				score += scoreSets(previousSlideTags, previousSlideNbTags, photo.tagBits, photo.nbTags)
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
					// Compute the score
					score += scoreSets(previousSlideTags, previousSlideNbTags, verticalSlideTags, verticalSlideNbTags)
					// Keep the current slide as the previous one
					previousSlideTags = verticalSlideTags
					previousSlideNbTags = verticalSlideNbTags
				}
			}
		}
	}

	return score
}

// scoreSets computes the score of two slides, each slide being a set of tags
func scoreSets(bits1 []uint64, nbBitsOn1 int, bits2 []uint64, nbBitsOn2 int) int {

	commonTags := 0
	for i := 0; i < len(bits1); i++ {
		tmp := bits1[i] & bits2[i]
		commonTags += bits.OnesCount64(tmp)
	}

	minimum := commonTags

	if nbBitSpecific1 := nbBitsOn1 - commonTags; nbBitSpecific1 < minimum {
		minimum = nbBitSpecific1
	}

	if nbBitSpecific2 := nbBitsOn2 - minimum; nbBitSpecific2 < minimum {
		minimum = nbBitSpecific2
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

	var breeder1, breeder2 int
	if useRoulette {
		breeder1 = selectBreederByRoulette(scores, -1)
		breeder2 = selectBreederByRoulette(scores, breeder1)
	} else {
		breeder1 = selectBreederByTournament(scores, -1)
		breeder2 = selectBreederByTournament(scores, breeder1)
	}
	return crossOver(population[breeder1], population[breeder2])
}

// selectBreederByRoulette selects randomly a breeder by using a roulette algorithm.
func selectBreederByRoulette(scores []int, forbiddenBreeder int) int {

	sumScores := 0
	for _, v := range scores {
		sumScores += v
	}

	// Rand limits are inclusive
	random := rand.Intn(sumScores)
	sum := 0

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
func selectBreederByTournament(scores []int, forbiddenBreeder int) int {

	tournamentSize := len(scores) / 8

	bestScore := -1
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
func readFile(fileName string) *[]photoDefinition {

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

		for j := 0; j < nbTags; j++ {
			tags[lineItems[j+2]] = true
		}

		photos[i] = photoDefinition{
			index:       i,
			orientation: photoOrientation,
			tags:        tags,
			tagBits:     nil,
			nbTags:      len(tags),
		}
	}

	// Collect all tags in a set
	tagsCollector := make(map[string]bool)
	// Collect all tags
	for _, photo := range photos {
		for tag := range photo.tags {
			tagsCollector[tag] = true
		}
	}

	// Get the key of the set to a list (ensure constant order)
	tags := make([]string, len(tagsCollector))
	idx := 0
	for tag := range tagsCollector {
		tags[idx] = tag
		idx++
	}

	// Number of uint64 needed by image
	nbUint64 := int(math.Ceil(float64(len(tags)) / 64.0))

	// Add the bytes to the image
	for i := 0; i < len(photos); i++ {
		photos[i].tagBits = make([]uint64, nbUint64)
	}

	// For each tag
	for tagIndex, tag := range tags {

		// Compute the position of the bit
		uint64idx := int(math.Floor(float64(tagIndex) / 64.0))
		posInUint64 := uint(tagIndex % 64)
		bitToSet := uint64(1) << posInUint64

		// For each photo
		for i := 0; i < len(photos); i++ {
			// If the photo has the tag
			if _, ok := photos[i].tags[tag]; ok {
				// Set the bit
				photos[i].tagBits[uint64idx] |= bitToSet
			}
		}
	}

	return &photos
}

// writeFile writes a genome to an output file
func writeFile(photos *[]photoDefinition, genome *[]int, fileName string) {

	file, err := os.Create(fileName)
	if err != nil {
		log.Panic(err)
	}

	defer func() {
		if e := file.Close(); e != nil {
			log.Panic(err)
		}
	}()

	previousVertical := -1

	for _, photoIndex := range *genome {

		photo := (*photos)[photoIndex]

		if photo.orientation == horizontal {
			_, err = file.WriteString(fmt.Sprintf("%v\n", photoIndex))
			if err != nil {
				log.Panic(err)
			}
		} else {
			if previousVertical == -1 {
				previousVertical = photoIndex
			} else {
				_, err = file.WriteString(fmt.Sprintf("%v %v\n", previousVertical, photoIndex))
				if err != nil {
					log.Panic(err)
				}
				previousVertical = -1
			}
		}
	}
}

// -------------------------------------------------------------------------------------------------------------------
//
//                                           UTILITIES
//
// -------------------------------------------------------------------------------------------------------------------

// analyzeScores displays information about the score computed (mean, avg, etc.)
func analyzeScores(scores []int) {

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
}
