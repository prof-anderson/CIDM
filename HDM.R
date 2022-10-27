#####################################################################################
### Functions : HDM, pairwise.comparison, pwc.matrix
### Script : HDM.R
### Description : Functions for determining weights from pairwise comparisons
### Author : Kevin Van Blommestein, 2014
#####################################################################################

#####################################################################################
### Setting up Environment
#####################################################################################
library(reshape2)
library(gtools)

#####################################################################################
### User defined functions
#####################################################################################
### Function to convert dataframe into square matrix required for
### pwc.weights function
pwc.matrix <- function(d, name.sep = "_", one.side = FALSE) {
	# Split column names by name.sep to determine criteria
	pw.names <- unlist(strsplit(colnames(d), split = name.sep))
	# Create data.frame with all pairwise comparisons
	df <- data.frame(c1 = pw.names[seq(1, length(pw.names), 2)],
									 c2 = pw.names[seq(2, length(pw.names), 2)], stringsAsFactors = FALSE)

	if (one.side)
		df <- rbind(df, data.frame(c1 = df[, 2], c2 = df[, 1]))

	#Check row names were specified correctly
	if ((nrow(df) != ncol(d) & !one.side) | (nrow(df) != 2 * ncol(d) & one.side))
		stop("Make sure column names are specified correctly. Use the correct name.sep for your column names.")

	# Determine unique criteria names
	c.names <- unique(as.vector(pw.names))

	if (is.null(row.names(d)))
		row.names(d) <- paste0("R_", 1:nrow(d))

	m <- list()
	# Iterate for each expert response and cast into square matrix format
	for (i in 1:nrow(d)){
		resp <- rownames(d)[i]
		if (one.side)
			di <- cbind(df, value = c(as.numeric(d[i, ]), 100 - as.numeric(d[i, ])))
		else
			di <- cbind(df, value = as.numeric(d[i, ]))

		m.i <- acast(di, c1 ~ c2, value.var="value", fill = 100)
		diag(m.i) <- 100
		# Append matrix to list of matrices for all expert responses
		m[[resp]] <- m.i
	}
	return(m)
}

### Function for calculating weights for pairwise comparisons
pwc.weights <- function(m){
	# Ensure data is in the correct format
	if (length(m) == 0)
		stop("List of matrices required, length(m) == 0")
	if (nrow(m[[1]]) != ncol(m[[1]]))
		stop("Ensure m is a list of square matrices")

	# Determine number of criteria
	n.c <- nrow(m[[1]])
	r.names <- names(m)

	# Save all results to matrices and return in list
	weights <- matrix(nrow = length(m), ncol = n.c, dimnames = list(names(m), row.names(m[[1]])))
	sdd <- matrix(nrow = length(m), ncol = n.c, dimnames = list(names(m), row.names(m[[1]])))
	inc <- matrix(nrow = length(m), ncol = 1, dimnames = list(names(m)))

	# Iterate for each expert response
	for (i in 1:length(m)){
		m.i <- m[[i]]
		# Determine number of criteria and criteria names
		n.c <- nrow(m.i)
		c.names <- colnames(m.i)

		# calculate ratio of parwise comparisons (pc11 / pc12, pc21 / pc22, etc.) - Matrix B
		m.b <- m.i / t(m.i)

		# Divide each column of matrix by other columns and find the column means - Matrix C
		p1 <- permutations(n.c, 2, 1:n.c)
		m.c <- colMeans(m.b[, p1[, 1]] / m.b[, p1[, 2]])

		# Determine each combination of multiplying m.c - Matrix D
		p2 <- permutations(n.c, n.c, v = 1:n.c)
		m.d <- data.frame(matrix(nrow = nrow(p2), ncol = ncol(p2)))
		colnames(m.d) <- colnames(p2) <- c.names

		# Calculate weights for each combination of m.c
		for (k in 1:nrow(p2)) {
			# Set last column criterion as 1 in m.d
			c.current <- p2[k, n.c]
			m.d[k, c.current] <- 1
			# Loop through remaining criteria and multiply with ratios from m.c
			# e.g (D * C/D * B/C * A/B), where D = 1 and ratios from m.c
			for (n in (n.c - 1):1) {
				c.next <- p2[k, n]
				m.d[k, c.next] <- m.d[k, c.current] * m.c[p1[, 1] == c.next & p1[, 2] == c.current]
				c.current <- c.next
			}
			# Normalise weights
			m.d[k, ] <- m.d[k, ]/sum(m.d[k, ])
		}

		# Calculate weights, standard deviations, and inconsistency
		weights[i, ] <- colMeans(m.d)
		sdd[i, ] <- apply(m.d, 2, sd) * sqrt((nrow(m.d) - 1) / nrow(m.d))
		inc[i, ] <- sqrt((1 / ncol(m.d)) * sum(sdd[i, ]^2))
	}

	weights.summ <- apply(weights, 2, function(x) c(mean(x), min(x), max(x), sd(x) * sqrt((length(x) - 1) / length(x))))
	row.names(weights.summ) <- c("mean", "min", "max", "sd")

	# Calculate disagreement among responses by first determining sd of weights, and then sd of sd's
	dis <- sd(apply(weights, 2, sd)  * sqrt((nrow(weights) - 1)/nrow(weights)))

	return(list(summary = list(weights = weights.summ, inconsistency = inc, disagreement = dis),
							detail = list(weights = weights, sd = sdd)))
}

### Convert data.frame to list for HDM function (currently only works for 2 levels)
df.list <- function(d, level.sep = ".", name.sep = "_", one.side = FALSE, levels = 2) {
	# Split column names into level identifier and pairwise identifier
	names.split <- strsplit(names(d), split = level.sep, fixed = TRUE)

	# Find parent names
	p.names <- sapply(names.split, "[[", 1)

	# Find which column is the objective
	p.len <- sapply(names.split, function(x) length(x))
	obj <- which(p.len == 1)

	# Make sure there is only one pairwise comparison at the objective level
	if (length(obj) < 1)
		stop("One objective required. Make sure all column names except the objective level have the level.sep character")

	# Find objective (highest) level criteria names
	obj.names <- unique(unlist(strsplit(p.names[obj], split = name.sep, fixed = TRUE)))

	# Create list of data to return
	lst <- vector("list", levels)

	# Add objective (highest) level to list
	lst[[1]] <- d[obj]

	d.sub <- d[, -obj]
	# Find all pairwise names
	pw.names <- sapply(names.split[-obj], "[[", 2)

	# Split each pairwise name at the highest level
	for (i in 1:length(obj.names)){
		# Select all columns from data for obj i
		lvl.name <- obj.names[i]
		lvl.cols <- which(p.names[-obj] %in% lvl.name)

		d.lvl <- d.sub[, lvl.cols]
		names(d.lvl) <- pw.names[lvl.cols]
		lst[[2]][[lvl.name]] <- d.lvl
	}

	return(lst)
}

### Calculate overall hdm weights
hdm <- function(d, name.sep = "_", levels = NA, one.side = FALSE) {
	# Determine number of levels in hierarchy
	if (is.data.frame(d)){
		d <- df.list(d, name.sep = name.sep, one.side = one.side)
		if (is.na(levels))
			stop("Specify number of levels in hierarchy")
	} else
		levels <- length(d)

	# There needs to be at least on level
	if (levels < 1)
		stop("At least one level in hierarchy is required")

	# Determine weights for main level (level 1)
	d.lvl1 <- d[[1]]
	m <- pwc.matrix(d.lvl1, name.sep = name.sep, one.side = one.side)
	pc <- pwc.weights(m)

	# Create list of results
	# r <- vector("list", levels)
	r <- d
	r[[1]] <- pc$summary

	if (levels > 1) {
		# Loop through each level
		for (lvl in 2:levels) {
			main.weights <- r[[lvl - 1]]$weights[1, ]

			r.lvl <- list()
			d.lvl <- d[[lvl]]

			# Check whether there is data for this level
			sub.levels <- length(d.lvl)
			if (sub.levels == 0)
				stop(paste("A least one pairwise comparison required for level", lvl))

			# Loop through each pairwise comparison for level lvl and save results
			for (sub.lvl in 1:sub.levels) {
				# Find weight for level above to calculate overall weights
				main.weight <- main.weights[match(names(d.lvl)[sub.lvl], names(main.weights))]
				m <- pwc.matrix(d.lvl[[sub.lvl]], name.sep = name.sep, one.side = one.side)
				pc <- pwc.weights(m)
				r[[lvl]][[sub.lvl]] <- pc$summary
				r[[lvl]][[sub.lvl]]$weights <- rbind(pc$summary$weights, obj = main.weight * pc$summary$weights[1,])
				r[[lvl]][[sub.lvl]]$disagreement <- sd(apply(main.weight * pc$detail$weights, 2, sd)  * sqrt((nrow(pc$detail$weights) - 1)/nrow(pc$detail$weights)))
			}
		}
	}
	return(r)
}


#####################################################################################
### Test
#####################################################################################

#-----------------------------------------------------------------------------------#
### Tests 1: Specify one-side of the pairwise comparisons, calculate weights + inconsistencies
#-----------------------------------------------------------------------------------#
# d1 consists of 5 expert responses, each row is an expert
d1 <- data.frame(A_B = c(90, 20, 60, 50, 50), A_C = c(10, 60, 80, 60, 50), B_C = c(90, 10, 30, 70, 50))
# Name rows with unique expert names
row.names(d1) <- paste0("R_", 1:5)

# Convert to square matrix required for pairwise comparison function
m <- pwc.matrix(d1, one.side = T)

pc <- pwc.weights(m)
# Expert 2 weights
pc


#-----------------------------------------------------------------------------------#
### Tests 2: Specify both sides of the pairwise comparisons, calculate weights + inconsistencies
#-----------------------------------------------------------------------------------#
# d1 consists of 5 expert responses, each row is an expert
d1 <- data.frame(A_B = c(90, 20, 60, 50, 50), B_A = c(10, 80, 40, 50, 50), A_C = c(10, 60, 80, 60, 50),
								 C_A = c(90, 40, 20, 40, 50), B_C = c(90, 10, 30, 70, 50), C_B = c(10, 90, 70, 30, 50))
# Name rows with unique expert names
row.names(d1) <- paste0("R_", 1:5)

# Convert to square matrix required for pairwise comparison function
m <- pwc.matrix(d1)

# Calculate weights
pc <- pwc.weights(m)
# Expert 2 weights
pc


#-----------------------------------------------------------------------------------#
### Tests 3: Multi-level pairwise comparisons (HDM) - Specify as List
#-----------------------------------------------------------------------------------#
# Create HDM model
# Model       Criteria  Subcriteria
#             -> A0   ->  A, B, C
# Objective -|
#             -> B0   ->  D, E, F

resp <- paste0("R_", 1:5)

level1 <- data.frame(A0_B0 = c(90, 40, 30, 40, 50), row.names = resp)
level2 <- list(A0 = data.frame(A_B = c(90, 20, 60, 50, 50), A_C = c(10, 60, 80, 60, 50), B_C = c(90, 10, 30, 70, 50), row.names = resp),
							 B0 = data.frame(D_E = c(90, 40, 30, 40, 50), D_F = c(10, 90, 70, 60, 50), E_F = c(10, 60, 70, 60, 50), row.names = resp))

d2 <- list(level1, level2)

r <- hdm(d2, one.side = TRUE)
r

#-----------------------------------------------------------------------------------#
### Tests 4: Multi-level pairwise comparisons (HDM) - Specify as dataframe
#-----------------------------------------------------------------------------------#

d <- data.frame(A_B = c(90, 20, 60, 50, 50), A_C = c(90, 20, 60, 50, 50), B_C = c(90, 20, 60, 50, 50),
								A.D_E = c(90, 40, 30, 40, 50), A.D_F = c(10, 90, 70, 60, 50), A.E_F = c(10, 60, 70, 60, 50),
								B.G_H = c(90, 40, 30, 40, 50), B.G_I = c(10, 90, 70, 60, 50), B.H_I = c(10, 60, 70, 60, 50),
								C.J_K = c(90, 40, 30, 40, 50), C.J_L = c(10, 90, 70, 60, 50), C.K_L = c(10, 60, 70, 60, 50))

dl <- df.list(d)

r <- hdm(dl, one.side = TRUE)

r


