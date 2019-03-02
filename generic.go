package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
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
	processFile("data/c_memorable_moments") // Mixed small
	//processFile("data/d_pet_pictures") // Mixed large
	//processFile("data/e_shiny_selfies") // All vertical

}

// populationSize is the size of the population.
const populationSize = 1000

// elitismFactor is the percentage of best parents that are kept in the next generation.
const elitismFactor = 0.1

// useRoulette defines if we use the roulette algorithm for breeding selection. Otherwise use tournament.
const useRoulette = true

// mutationFactor is the factor applied to each gene to determine if it mutates, keep it small
const mutationFactor = 0.001

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
}

// scoreResult glues together an index in the population and its score.
type scoreResult struct {
	genomeIndex int
	score       int
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

	messages := make(chan scoreResult)

	population := createInitialPopulation(len(*photos), populationSize)
	scores := make([]int, populationSize)

	for round := 0; round < 100; round++ {

		fmt.Printf("Generation: %v\n", round)

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

		fmt.Printf("genome scored in %v\n", time.Since(start))

		analyzeScores(scores)

		start = time.Now()
		population = breed(population, scores)

		fmt.Printf("population breeded in %v\n", time.Since(start))

		start = time.Now()
		mutate(population)

		fmt.Printf("population mutated in %v\n", time.Since(start))
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

// scoreGenome scores a genome. The resulting score "should" be the Google Score. The slideware is generated as a slide
// for each horizontal photo and a slice for each pair of vertical photo. For pair the slide is generated when the
// second member of the pair is found
func scoreGenome(photos *[]photoDefinition, genome *[]int) int {

	var previousVertical *photoDefinition
	var previousSliceTags map[string]bool

	score := 0

	for _, photoIndex := range *genome {

		photo := (*photos)[photoIndex]

		// If initialization of scoring
		if previousSliceTags == nil {
			if photo.orientation == horizontal {
				previousSliceTags = photo.tags
			} else {
				if previousVertical == nil {
					previousVertical = &photo
				} else {
					previousSliceTags = addSets(previousVertical.tags, photo.tags)
					previousVertical = nil
				}
			}
		} else {
			if photo.orientation == horizontal {
				score += scoreSets(previousSliceTags, photo.tags)
				previousSliceTags = photo.tags
			} else {
				if previousVertical == nil {
					previousVertical = &photo
				} else {
					tags := addSets(previousVertical.tags, photo.tags)
					score += scoreSets(previousSliceTags, photo.tags)
					previousSliceTags = tags
				}
			}
		}
	}

	return score
}

// addSets adds two sets into a new one. This function is used to for having the tags of a pair of vertical photos.
func addSets(set1 map[string]bool, set2 map[string]bool) map[string]bool {
	results := make(map[string]bool)
	for k := range set1 {
		results[k] = true
	}
	for k := range set2 {
		results[k] = true
	}
	return results
}

// scoreSets computes the score of two slides, each slide being a set of tags
func scoreSets(set1 map[string]bool, set2 map[string]bool) int {

	commonTags := 0
	for tag := range set1 {
		if _, ok := set2[tag]; ok {
			commonTags++
		}
	}

	minimum := commonTags

	if v := len(set1) - commonTags; v < minimum {
		minimum = v
	}

	if v := len(set2) - commonTags; v < minimum {
		minimum = v
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

	// Create new children
	for i := elitismKept; i < len(population); i++ {

		var breeder1, breeder2 int
		if useRoulette {
			breeder1 = selectBreederByRoulette(scores)
			breeder2 = selectBreederByRoulette(scores)
		} else {
			breeder1 = selectBreederByTournament(scores)
			breeder2 = selectBreederByTournament(scores)
		}
		children := crossOver(population[breeder1], population[breeder2])
		results[i] = children
	}

	return results
}

// selectBreederByRoulette selects randomly a breeder by using a roulette algorithm.
func selectBreederByRoulette(scores []int) int {

	sumScores := 0
	for _, v := range scores {
		sumScores += v
	}

	// Rand limits are inclusive
	random := rand.Intn(sumScores)
	sum := 0
	for i, v := range scores {
		sum += v
		if sum >= random {
			return i
		}
	}
	// should not happen
	return 0
}

// selectBreederByRoulette selects randomly a breeder by using a tournament algorithm.
func selectBreederByTournament(scores []int) int {

	tournamentSize := len(scores) / 8

	bestScore := -1
	var bestBreeder int

	for i := 0; i < tournamentSize; i++ {
		random := rand.Intn(len(scores))
		if scores[random] > bestScore {
			bestScore = scores[random]
			bestBreeder = random
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

	for i := start; i < length; i++ {
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
//                                    MUTATTION FUNCTIONS
//
// -------------------------------------------------------------------------------------------------------------------

// mutate applies random mutation on a population. Mutations are made by swapping the mutated gene with another.
func mutate(population []*[]int) {

	// For all genomes
	for _, genome := range population {
		genomeSize := len(*genome)
		// For all genes
		for i, gene := range *genome {
			// If the gene is to be mutated
			if rand.Float64() < mutationFactor {

				// Make a random swap
				idxTarget := rand.Intn(genomeSize)
				temp := (*genome)[idxTarget]
				(*genome)[idxTarget] = gene
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

		nbTags, err := strconv.Atoi(lineItems[1])
		if err != nil {
			log.Fatal(err)
		}

		tags := make(map[string]bool, nbTags)

		for j := 0; j < nbTags; j++ {
			tags[lineItems[j+2]] = true
		}

		photos[i] = photoDefinition{
			index:       i,
			orientation: photoOrientation,
			tags:        tags,
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
