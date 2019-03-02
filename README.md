#Google Hash Code 2019 - Genetic algorithm implementation

This project tries to solve the Google [Hash Code 2019](https://codingcompetitions.withgoogle.com/hashcode/) problem 
using a genetic algorithm for searching the best solution.

Status: This approach is not working well at all... I just keep it for the sake of memory.

There are two main issues:

* The number of combination is very high (up to *factorial(80.000)*), so the initial population for searching should be 
also large. Unfortunately, the time for scoring each solution may be quite high (currently ~ 80 ms on my computer), so 
a step of the algorithm takes a lot of time
* Although the number of combination is quite large, the difference between two combinations may be quite small. So 
finding an optimal solution requires a lot of steps ... that are time consuming. 

As of now, just the scoring function is run in parallel as it is the most time consuming operation

#Versions

 * v0.0.1: Original version
 * v0.0.2: Enhance peformances (up to 30x)

# License

Copyright 2019 Thomas Wuillemin  <thomas.wuillemin@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this project or its content except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.