#!/usr/bin/env Rscript

.tinyOver = function(
	x,                # Pathway ID
	p.counts.total,   # Table of genes for all pathways
	p.counts.sig,     # Table of significant genes for all pathways
	g.sig,            # Number of significant genes
	g.total,          # Number of genes
	stat,             # Statistic
	alt               # Alternative hypothesis
) {

	x = as.character(x)
	
	x.total = as.integer(p.counts.total[x])  # Number of genes in pathway
	x.sig = as.integer(p.counts.sig[x])      # Number of genes in pathway AND foreground
	
	if (stat %in% c('fisher','chisq')) {
		contingency.matrix = matrix(c(
			x.sig,                               # In pathway AND foreground
			x.total - x.sig,                     # In pathway but not foreground
			g.sig - x.sig,                       # In foreground but not pathway
			g.total - x.total - g.sig + x.sig    # In neither foreground nor pathway
		),nrow=2,dimnames=list(                  # Dimension names are not necessary,
			c('A+B+','A+B–'),                    #    but are useful for debugging
			c('A–B+','A–B–')                     #  
		))            
		
		test.result = if (stat == 'fisher') {
			fisher.test(contingency.matrix,alternative=alt)
		} else if (stat == 'chisq') {
			suppressWarnings(chisq.test(contingency.matrix))
		}
	} else if (stat == 'binom') {
		test.result = binom.test(x.sig,x.total,p = g.sig/g.total,alternative=alt)
	}
	if (stat == 'fisher') {
		out = data.frame(
			pathway_id = x,
			estimate = test.result$estimate,
			pval = test.result$p.value,
			annotated = x.total,
			significant = x.sig,
			expected = x.sig / test.result$estimate
		)
	} else if (stat == 'binom') {
		out = data.frame(
			pathway_id = x,
			statistic = test.result$statistic,
			estimate = test.result$estimate,
			pval = test.result$p.value,
			annotated = x.total,
			significant = x.sig,
			expected = x.total * g.sig/g.total
		)
	} else if (stat == 'chisq') {
		out = data.frame(
			pathway_id = x,
			statistic = test.result$statistic,
			pval = test.result$p.value,
			annotated = x.total,
			significant = x.sig,
			expected = test.result$expected[1]
		)
	}
	out
}

.tinyEnrich = function(
	x,                # Pathway ID
	p2g,              # List of genes assigned to pathways
	v,                # Test vector
	stat,             # Statistic
	alt               # Alternative hypothesis
) {

	x = as.character(x)
	
	in.genes = p2g[[x]]
	x.total = length(in.genes)               # Number of genes in pathway

	in.pathway = v[in.genes]
	ex.pathway = v[setdiff(names(v),in.genes)]

	test.result = if (stat == 'ks') {
		ks.test(in.pathway,ex.pathway,alternative=alt)
	} else if (stat == 'wilcox') {
		wilcox.test(in.pathway,ex.pathway,alternative=alt)
	} else if (stat == 't') {
		t.test(in.pathway,ex.pathway,alternative=alt)
	}
	
	if (stat == 'ks') {
		out = data.frame(
			pathway_id = x,
			statistic = test.result$statistic,
			pval = test.result$p.value,
			annotated = x.total
		)
	} else if (stat == 'wilcox') {
		out = data.frame(
			pathway_id = x,
			statistic = test.result$statistic,
			pval = test.result$p.value,
			annotated = x.total
		)
	} else if (stat == 't') {
		out = data.frame(
			pathway_id = x,
			statistic = test.result$statistic,
			estimate = diff(test.result$estimate),
			pval = test.result$p.value,
			annotated = x.total
		)
	}
	out
}

#' Simple (*prim*itive and *prim*ate-inspired) enrichment analysis workflow
#'
#' This function provides a workflow for running simple enrichment tests using
#' annotations and gene lists or gene scores.
#' 
#' @param test.values A vector of test values. May be a character, logical, integer, or numeric class. If a character vector is given, the `background` argument is required.
#' @param annotation Either an object of "pathwayAnnotation" class or a named list where names are gene IDs and contents are pathway IDs. If a list is given, the `definition` argument is required.
#' @param alternative A character string specifying the alternative hypothesis, must be one of `"greater"` (default), `"less"`, or `"two.sided"`.
#' @param background A character vector of background genes. Used only when `test.values` is a character vector, ignored otherwise.
#' @param statistic Character string specifing which test to use. Possible values are `"fisher"`, `"binom"`, `"chisq"`, `"ks"`, `"wilcox"`, and `"t"`, 
#' @param node.size A number representing the minimum number of genes per pathway. The default is `10`.
#' @param definition A data frame with the columns `pathway_id` and `pathway_name`, along with additional optional columns. Ignored if `annotation` is of the pathwayAnnotation class.
#' @param n.cores A number specifying the number of cores to use.
#' @return A data frame of pathway IDs, test statistics, and associated pathway metadata.
#' @export
primEnrich = function(    # Simple enrichment test
	test.values,          # character vector of differentially expressed genes OR
	                      #     a named vector of class integer, numeric, or logical
	annotation,           # Either object of "pathwayAnnotation" class or a named list
	                      #     where names are genes and contents are pathway IDs.
	alternative='greater',# Direction for test ("greater" indicates alternative hypothesis
	                      #     that pathway-associated values are more positive)
	background=NULL,      # character vector of background genes (all possible genes);
	                      #     non-DE genes is also ok as the background would be the
	                      #     union of the two.
	statistic='fisher',   # type of test to run. Possible values: 'fisher','chisq',
	                      #     'binom','ks','wilcox','t'
	node.size=10,         # nodeSize argument for topGO, prunes the hierarchy of terms
	                      #     that have fewer than this number of genes
	definition=NULL,      # data frame consisting of columns pathway_id and pathway_name
	                      #     Additional columns may be optionally included. Ignored
	                      #     if annotation is of class "pathwayAnnotation". If NULL,
	                      #     pathway names are left blank.
	n.cores = NULL
) {
	if (statistic %in% c('fisher','chisq','binom')) {
		test.type = 'overrepresentation' 
	} else if (statistic %in% c('ks','wilcox','t')) {
		test.type = 'enrichment'
	} else {
		stop(paste0('Unrecognized statistic: ',statistic))
	}
	
	if (statistic == 'chisq') {
		if (!is.null(alternative)) warning('alternative is ignored for chi-squared tests.')
		alternative = '' # set to a dummy value to avoid NULL issues
	} else if (is.null(alternative)) {
		stop(paste0('alternative cannot be NULL when statistic = "',statistic,'"'))
	}
	if (statistic %in% c('fisher','binom') & alternative != 'greater') {
		warning(paste0('Running "',statistic,'" test with alternative = "',alternative,'" rather than alternative = "greater". This is an atypical usage of this test.'))
	}
	
	if (test.type == 'overrepresentation') {
		if (length(intersect(class(test.values),c('character')))) {
			foreground = test.values
			if (is.null(background)) stop('background may not be NULL when test.values is of class "character"')
			if (!'character' %in% class(background)) stop('background must be of class "character"')
			if (!all(foreground %in% background)) {
				message('Some test genes were not in background and are being added.')
				background = union(foreground,background)
			}
			test.vector = background %in% foreground
			names(test.vector) = background
		} else {
			if (is.null(names(test.values))) stop('test.values must be a named vector.')
			if (length(intersect(class(test.values),c('numeric','integer','logical')))) {
				test.vector = as.logical(test.values)
				names(test.vector) = names(test.values)
			} else {
				stop('test.values is not of a recognized class')
			}
		}
	} else {
		if ('character' %in% class(test.values)) stop(paste0('test.values cannot be of class "character" when statistic = ',statistic))
		if (is.null(names(test.values))) stop('test.values must be a named vector.')
		test.vector = as.numeric(test.values)
		names(test.vector) = names(test.values)
	}
	
	# if (!is.null(n.cores)) suppressPackageStartupMessages(require(parallel))
	
	if (class(annotation) == 'pathwayAnnotation') {
		gene2pwy = annotation@annotation
		definition = annotation@definition
	} else {
		gene2pwy = annotation
		if (is.null(definition)) {
			warning('No definition object provided. Pathways will be named by pathway IDs, which may or may not be interpretable to human eyes.')
			pathway.ids = as.character(sort(unique(unlist(annotation))))
			definition = data.frame(
				pathway_id = pathway.ids,
				pathway_name = pathway.ids
			)
		} else if (!'data.frame' %in% class(definition)) {
			stop('Argument definition must be a data frame.')
		}
	}
	if (sum(duplicated(definition$pathway_id))) {
		stop('pathway_id in definition data frame cannot contain duplicates.')
	}
	rownames(definition) = definition$pathway_id
	
	# Keep only annotations of genes in data
	gene2pwy.genes = gene2pwy[intersect(names(test.vector),names(gene2pwy))]
	
	gene2pwy.df = do.call(rbind,if (is.null(n.cores)) {
		lapply(
			names(gene2pwy.genes),
			function(x)	data.frame(gene_id = x,pathway_id = sort(unique(gene2pwy.genes[[x]])))
		)
	} else {
		parallel::mclapply(
			names(gene2pwy.genes),
			function(x)	data.frame(gene_id = x,pathway_id = sort(unique(gene2pwy.genes[[x]]))),
			mc.cores = n.cores
		)
	})
	
	# Node filter
	gene2pwy.df = subset(gene2pwy.df,pathway_id %in% names(which(table(gene2pwy.df$pathway_id) > node.size)))
	
	# Keep only test genes that overlap
	test.vector = test.vector[intersect(names(test.vector),gene2pwy.df$gene_id)]
	
	# Add test values to data frame
	gene2pwy.df$test = test.vector[gene2pwy.df$gene_id]
	
	pathway.list = sort(unique(as.character(gene2pwy.df$pathway_id)))
	
	if (test.type == 'overrepresentation') {
		pwy.counts = table(gene2pwy.df$pathway_id)     # Total annotated (per pathway)
		pwy.foreground.counts = table(factor(          # Total significant (per pathway)
			subset(gene2pwy.df,test)$pathway_id,       #
			levels=names(pwy.counts)                   #
		))                                             #
		n.foreground = sum(test.vector)                # Total significant (overall)
		n.background = length(test.vector)             # Total genes
		
		result = do.call(rbind,if (is.null(n.cores)) {
			lapply(pathway.list, .tinyOver, p.counts.total = pwy.counts, p.counts.sig = pwy.foreground.counts, g.sig = n.foreground, g.total = n.background, stat = statistic, alt = alternative)
		} else {
			parallel::mclapply(pathway.list, .tinyOver, p.counts.total = pwy.counts, p.counts.sig = pwy.foreground.counts, g.sig = n.foreground, g.total = n.background, stat = statistic, alt = alternative, mc.cores = n.cores)
		})
	} else if (test.type == 'enrichment') {
		pwy2gene = split(gene2pwy.df$gene_id,gene2pwy.df$pathway_id)
		
		if (statistic == 'ks' & alternative %in% c('less','greater')) {
			warning('Note that test is being run with an intuitive "alternative" argument, but one that differs from the argument to ks.test. See documentation for more details.')
			if (alternative == 'greater') {
				alternative = 'less'
			} else {
				alternative = 'greater'
			}
		}
		
		result = do.call(rbind,if (is.null(n.cores)) {
			lapply(pathway.list, .tinyEnrich, p2g = pwy2gene, v = test.vector, stat = statistic, alt = alternative)
		} else {
			parallel::mclapply(pathway.list, .tinyEnrich, p2g = pwy2gene, v = test.vector, stat = statistic, alt = alternative, mc.cores = n.cores)
		})
	}
	
	result = merge(result,definition,by='pathway_id',all.x=TRUE)
	
	result$padj = p.adjust(result$pval,'fdr')
	result = result[order(result$pval),]
	rownames(result) = NULL
	
	# Return data frame in order
	required.columns = intersect(c('pathway_id','pathway_name','annotated','significant','expected','estimate','statistic','pval','padj'),colnames(result))
	drop.columns = c('GO.ID','Term')
	result[,c(required.columns,setdiff(colnames(result),c(required.columns,drop.columns)))]
}
