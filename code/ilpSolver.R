ilpSolver<-function(g, K, SETS)
{
	# To run the example
	#g <- graph( c(1,2, 4,3, 2,3, 1,3, 4,1 ,2,4, 1,5, 2,5, 3,5, 4,5 ), n=4 ) # complete graph
  #E(g)$weight=0.8 # weights for edges - should be Ws
	#E(g)[3]$weight=0.5
	#E(g)[6]$weight=1
	#E(g)[2]$weight=0.3
	#E(g)[7]$weight=0.9
	#E(g)[9]$weight=0
	#E(g)[10]$weight=0.4

	#input:  complete weighted graph, K = size of group, SETS = number of sets to retreive
	#output: matrix of chosen sets, ordered from the set the scored the best to the least
	#SETS = 15
	#K = 3

	V = gorder(g)
	E = ecount(g)

	# Maximumn manners to choose K vertices out of the g
	SETS_maximum = min(choose(V, K), choose(E,(K-1)*K/2))

	# Take the minimum of either maximum possible sets to choose from, or the number of sets requested
	SETS = min(SETS_maximum,SETS)
	chosen_SETS = matrix(nrow = SETS, ncol = V)
	score_SETS = rep(-1,SETS)
	status_SETS = rep(-1,SETS)
	
	# Should be done only once:

	# Open cplex enviorement
	env <- openEnvCPLEX()
	# Create a new problem
	prob <- initProbCPLEX(env)

	# Set the  problem:

	# Number of variables
	nc = V + E
	# Number of constraints (6,7,8 are per edge, and one more to limit the number of edges and to limit the number of edges + one for each set)
	nr = 3*E + 2 + SETS

	# Number of non zero elements: constraint 6 has edge + vertex, constraint 7 has edge + vertex, constraint 7 has edge + 2 vertices, and all three constraints are repeated for each edge
	# constraint 5 has V elements and constraint 10 has E elements
	nz = 7*E + V + E + V * SETS

	# Objective function: all chosen edges contribute their Ws scores
	obj <- c(rep(0,V), E(g)$weight)

	# Set right hand set: 6 + 7 garantie that if an edge is chosen in the Set then both its nodes are chosen. 8 garanties that if both nides are chosen then thhe edge is chosen
	# Right hand side, 5 that exactly K vertices are chosen and 10 that exactly  (K*(K-1))/2 edges are chosen, and finally for each set that the previous chosen nodes are not rechosen
	rhs <- c(rep(c(0,0,1),E),K, (K*(K-1))/2, rep(K-1, SETS))

	# Sense correlates to the sense (=, <=, or =>) of constraints: 6,7,8 in the ILP formula are <= and |chosen_V| == K, and |chosen_E| == K*(K-1)/2
	#sense <- c(rep("L",3*E),"E", "E", rep("L",SETS))
	sense <- c(rep("L",3*E),"E", "L", rep("L",SETS)) # try not to force a complete set

	# All variables are binary integers
	ctype <- rep("I", V+E)
	lb <- rep(0.0, V+E)
	ub <- rep(1.0, V+E)
	
	# Where each variable begin, vertices first then edges.
	# Vertices:  for each edge, 4 vertices participate in its contraints (if edge (u,v) choose u(1), if edge (u,v) choose v(2), if u(3) and v(4) are chosen choose (u,v))
	# + one for |chosen_V| == K + one for each set
	# and each edge has 3 constraints (for (u,v): if edge (u,v) choose u, if edge (u,v) choose v, if u and v are chosen choose (u,v)) + one for |chosen_E| == K*(K-1)/2
	
	beg1 <- c(rep(-1,V))
	cnt <- c(rep (0,V), rep (4, E))
	
	# each vertex participates in:
	# 2 constraints for each edge + one for |chosen_V| == K + one for each set
	prevBeg1 <- 0
	for (vertex in 1:V) {
	  dVertex = degree(g, vertex, mode="all")
	  
	  cnt[vertex] <- 2*dVertex + 1 +SETS
	  beg1[vertex] <- prevBeg1
	  prevBeg1 <- prevBeg1 + cnt[vertex]
	}
	
	# beginning for edges
	beg2 <- seq (prevBeg1,by = 4, length.out = E)
	beg <- c(beg1,beg2)
	
	# OLD and not clear why it ever worked:
	#constraints4vertex = 2*(V - 1) + 1 + SETS
	#beg1 <- seq (0,by = constraints4vertex, length.out = V)
	#beg2 <- seq (beg1[V] + constraints4vertex,by = 4, length.out = E)
	#beg <- c(beg1,beg2)
	#cnt <- c(rep (constraints4vertex, V), rep (4, E))

	# Indexes: which variable participates in which constraints
	ind <- rep(-1,7*E + V + E)
	val <- rep(0,7*E + V + E)
	edges = get.edges(g, E(g))
	i = 0
	beg_edges = beg[V+1] + 1
	cntCopy <- rep (0,V)

	# Constraints 6,7,8:
	# 6. if (u,v) is chosen, then choose u
	# 7. if (u,v) is chosen, then choose u
	# 8. if u and v are chosen, choose (u,v)
	for (edge in 1:E) {
	  
	  curr_edge = beg_edges + (edge - 1) * 4
	  u = edges[edge,1]
	  v = edges[edge,2]
	  
	  # if (u,v) is chosen, then choose u
	  ind[curr_edge] = i
	  ind[1 + beg[u] + cntCopy[u]] = i
	  val[curr_edge] = 1
	  val[1 + beg[u] + cntCopy[u]] = -1
	  i = i + 1
	  cntCopy[u] = cntCopy[u] + 1
	  
	  # if (u,v) is chosen, then choose u
	  ind[curr_edge + 1] = i
	  ind[1 + beg[v] + cntCopy[v]] = i
	  val[curr_edge + 1] = 1
	  val[1 + beg[v] + cntCopy[v]] = -1
	  i = i + 1
	  cntCopy[v] = cntCopy[v] + 1
	  
	  # if u and v are chosen, choose (u,v)
	  ind[curr_edge + 2] = i
	  ind[1 + beg[u] + cntCopy[u]] = i
	  ind[1 + beg[v] + cntCopy[v]] = i
	  val[curr_edge + 2] = -1
	  val[1 + beg[u] + cntCopy[u]] = 1
	  val[1 + beg[v] + cntCopy[v]] = 1
	  cntCopy[u] = cntCopy[u] + 1
	  cntCopy[v] = cntCopy[v] + 1
	  i = i + 1
	}

	# Constraint 5. Choose only K of the vertices
	for (vertex in 1:V) {
	  
	  ind[1 + beg[vertex] + cntCopy[vertex]] = i
	  val[1 + beg[vertex] + cntCopy[vertex]] = 1
	  cntCopy[vertex] = cntCopy[vertex] + 1
	}
	i = i + 1

	# Constraint 10. Choose only K(K-1)/2 of the vertices
	for (edge in 1:E) {
	  curr_edge = beg_edges + (edge - 1) * 4
	  ind[curr_edge + 3] = i
	  val[curr_edge + 3] = 1
	}
	i = i + 1
	
	# Constraint 9. Prepare SETS constraints: for each set the chosen vertices will be later set to 1 
	for (set in 1:SETS) {
	  for (vertex in 1:V) {
		
  		ind[1 + beg[vertex] + cntCopy[vertex]] = i
  		val[1 + beg[vertex] + cntCopy[vertex]] = 0
  		cntCopy[vertex] = cntCopy[vertex] + 1
	  }
	  i = i + 1
	}

	# Get solution from CPLEX
	copyLpCPLEX(env, prob, nc, nr, CPX_MAX, obj, rhs, sense,  beg, cnt, ind, val, lb, ub)
	copyColTypeCPLEX(env, prob, ctype)

	writeProbCPLEX(env, prob, "prob.lp")

	mipoptCPLEX(env, prob)
	solution = solutionCPLEX(env, prob)
	print(solution$lpstat)
	
	# if successful
	if (solution$lpstat == 101 | solution$lpstat == 102) {
		  chosen_SETS[1,] <- solution$x[1:V]
		  score_SETS[1] = solution$objval
		  status_SETS[1] = solution$lpstat
	}

	# Add a constraint limiting the number of vertices chosen
	SET_start <- 3*E
	curr_Solution <- 2
	while (curr_Solution <= SETS) {
		chgCoefListCPLEX(env, prob, V, rep(SET_start + curr_Solution,V), seq (0,by = 1, length.out = V), solution$x[1:V])
	  
	  # Find current set
	  # If cplex failed to find a new set it will throw an error (althohg this is a valis status in this case, as sets may be simply exhausted)
	    # Try to solve problem
	    mipoptCPLEX(env, prob)
	    writeProbCPLEX(env, prob, "prob.lp")
	    
	    solution <- try(solutionCPLEX(env, prob))

	    if (class(solution) == "list") { 

        # if successful
        print(solution$lpstat)
        if (solution$lpstat == 101 | solution$lpstat == 102) {

          chosen_SETS[curr_Solution,] <- round(solution$x[1:V]) # Cplex suppose to assign 1 / 0, but due to accuracy issues it may return values close to 1 or 0
          score_SETS[curr_Solution] = solution$objval
          status_SETS[curr_Solution] = solution$lpstat
          curr_Solution <- curr_Solution + 1
        }
	    }
	}
	
	# Free memory, allacated to the problem object
	delProbCPLEX(env, prob)
	closeEnvCPLEX(env)

	return (list("sets" = chosen_SETS, "scores" = score_SETS, "status" = status_SETS))
}
