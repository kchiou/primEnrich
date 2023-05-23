#!/usr/bin/env Rscript

#' Wrapper for Fisher's exact test implemented in topGO
#'
#' This function provides a workflow for running topGO enrichment tests which
#' take into account the topology of the GO graph.
#' 
#' @param test.values A vector of test values. May be a character, logical, integer, or numeric class. If a character vector is given, the `background` argument is required.
#' @param annotation Either an object of "pathwayAnnotation" class or a named list where names are gene IDs and contents are pathway IDs. If a list is given, the `definition` argument is required.
#' @param background A character vector of background genes. Used only when `test.values` is a character vector, ignored otherwise.
#' @param node.size A number representing the minimum number of genes per pathway. The default is `10`.
#' @param definition A data frame with the columns `pathway_id` and `pathway_name`, along with additional optional columns. Ignored if `annotation` is of the pathwayAnnotation class.
#' @param algorithm Character string specifing which algorithm to use. Possible values are given by `topGO::whichAlgorithms()`.
#' @param go.ontology A character string specifying the namespace for which to limit GO annotations. Must be one or more of `"BP"`, `"CC"`, or `"MF"`.
#' @return A data frame of GO pathway IDs, test statistics, and associated metadata.
#' @export
topgoFisher = function(   # Fisher's exact test with topGO
	test.values,          # character vector of differentially expressed genes
	annotation,           # Either object of "pathwayAnnotation" class or a named list
	                      #     where names are genes and contents are GO terms.
	background=NULL,      # character vector of background genes (all possible genes);
	                      #     non-DE genes is also ok as the background would be the
	                      #     union of the two.
	node.size=10,         # nodeSize argument for topGO, prunes the hierarchy of terms
	                      #     that have fewer than this number of genes
	definition=NULL,      # data frame consisting of columns pathway_id and pathway_name
	                      #     Additional columns may be optionally included. Ignored
	                      #     if annotation is of class "pathwayAnnotation". If NULL,
	                      #     pathway names are left blank.
	algorithm='weight01', # algorithm for Fisher's exact test, one of: c("classic",
	                      #     "elim", "weight", "weight01", "lea", "parentchild")
	go.ontology='BP'      # GO ontology, if restricting the enrichment to a particular
	                      #     domain.
) {
	
	if (class(test.values) == 'character') {
		foreground = test.values
		if (is.null(background)) stop('background may not be NULL when foreground is of class character')
		if (!all(foreground %in% background)) {
			warning('Some foreground genes were not in background and are being added')
			background = union(foreground,background)
		}
		test.vector = as.integer(background %in% foreground)
		names(test.vector) = background
	} else if (length(intersect(class(test.values),c('numeric','integer','logical')))) {
		test.vector = as.integer(as.logical(test.values))
		names(test.vector) = names(test.values)
		foreground = names(which(test.vector>0))
		background = names(test.vector)
	}
	
	suppressPackageStartupMessages(require(topGO))
	
	if (class(annotation) == 'pathwayAnnotation') {
		if (!grepl('^GO *',annotation@source)) warning('Slot "source" of pathwayAnnotation object should be "GO"')
		gene2go = annotation@annotation
		definition = annotation@definition
	} else {
		gene2go = annotation
	}
	
	test.object = suppressMessages(new('topGOdata',
		description = 'Simple session',
		ontology = go.ontology,
		allGenes = test.vector,
		geneSelectionFun = function(x) x > 0,
		nodeSize = node.size,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	))
	test.result = suppressMessages(topGO::runTest(test.object,statistic='fisher',algorithm=algorithm))
	result = topGO::GenTable(
		test.object,
		pval=test.result,
		topNodes=length(test.result@score),
		numChar=1000 # Set ridiculously high to ensure that the full name is printed
	)
	result$go_id = as.character(result$GO.ID)
	result$go_name = result$Term
	result$annotated = result$Annotated
	result$significant = result$Significant
	result$expected = result$Expected
	result$estimate = with(result,significant/expected)
	result$pval = test.result@score[result$go_id] # Rewrite p value to make numeric
	result$padj = p.adjust(result$pval,'fdr') # Calculate FDR
	
	# If the definition object exists has additional columns (besides pathway_id and pathway_name), tack them on
	# if (exists('definition',inherit=FALSE)) {
	if (!is.null(definition)) {
		definition.join = definition
		definition.join$go_id = definition.join$pathway_id
		definition.join$pathway_id = NULL
		definition.join$pathway_name = NULL
		result = merge(result,definition.join,by='go_id',all.x=TRUE)
	}
	
	result = result[order(result$pval),]
	rownames(result) = NULL
	
	required.columns = c('go_id','go_name','annotated','significant','expected','estimate','pval','padj')
	drop.columns = c('GO.ID','Term','Annotated','Significant','Expected')
	result[,c(required.columns,setdiff(colnames(result),c(required.columns,drop.columns)))]
}

#' Wrapper for two-sample Kolmogorov-Smirnov test implemented in topGO
#'
#' This function provides a workflow for running topGO enrichment tests, which
#' take into account the topology of the GO graph.
#' 
#' @param test.values A named vector of test values. Must be coercable to numeric class.
#' @param annotation Either an object of "pathwayAnnotation" class or a named list where names are gene IDs and contents are pathway IDs. If a list is given, the `definition` argument is required.
#' @param alternative A character string specifying the alternative hypothesis, must be either `"less"` (default) or `"greater"`.
#' @param node.size A number representing the minimum number of genes per pathway. The default is `10`.
#' @param definition A data frame with the columns `pathway_id` and `pathway_name`, along with additional optional columns. Ignored if `annotation` is of the pathwayAnnotation class.
#' @param algorithm Character string specifing which algorithm to use. Possible values are given by `topGO::whichAlgorithms()`.
#' @param go.ontology A character string specifying the namespace for which to limit GO annotations. Must be one or more of `"BP"`, `"CC"`, or `"MF"`.
#' @return A data frame of GO pathway IDs, test statistics, and associated metadata.
topgoKS = function(       # Fisher's exact test with topGO
	test.values,          # Named vector of integer or numeric class. Values may be, for
	                      #     example, p values or test statistics
	annotation,           # Either object of "pathwayAnnotation" class or a named list
	                      #     where names are genes and contents are pathway IDs.
	alternative='less',   # Direction for test ("greater" indicates alternative hypothesis
	                      #     that pathway associated values are more positive)
	node.size=10,         # nodeSize argument for topGO, prunes the hierarchy of terms
	                      #     that have fewer than this number of genes
	definition=NULL,      # data frame consisting of columns pathway_id and pathway_name
	                      #     Additional columns may be optionally included. Ignored
	                      #     if annotation is of class "pathwayAnnotation". If NULL,
	                      #     pathway names are left blank.
	algorithm='weight01', # algorithm for Fisher's exact test, one of: c("classic",
	                      #     "elim", "weight", "weight01", "lea", "parentchild")
	go.ontology='BP'      # GO ontology, if restricting the enrichment to a particular
	                      #     domain.
) {
	
	if (is.null(names(test.values))) {
		stop('Argument test.values must be named according to gene IDs.')
	}
	
	gene.names = names(test.values)
	if (!(('integer' %in% class(test.values)) || ('numeric' %in% class(test.values)))) {
		stop('Argument test.values must be numeric.')
	}
	
	test.vector = as.numeric(test.values)
	names(test.vector) = gene.names
	
	suppressPackageStartupMessages(require(topGO))
	
	if (class(annotation) == 'pathwayAnnotation') {
		if (!grepl('^GO *',annotation@source)) warning('Slot "source" of pathwayAnnotation object should be "GO"')
		gene2go = annotation@annotation
		definition = annotation@definition
	} else {
		gene2go = annotation
	}
	
	test.object = suppressMessages(new('topGOdata',
		description = 'Simple session',
		ontology = go.ontology,
		allGenes = if (alternative=='less') test.vector else -test.vector,
		geneSelectionFun = function(x) x > -Inf,
		nodeSize = node.size,
		annotationFun = annFUN.gene2GO,
		gene2GO = gene2go
	))
	test.result = suppressMessages(topGO::runTest(test.object,statistic='ks',algorithm=algorithm))
	result = topGO::GenTable(
		test.object,
		pval=test.result,
		topNodes=length(test.result@score),
		numChar=1000 # Set ridiculously high to ensure that the full name is printed
	)
	result$go_id = as.character(result$GO.ID)
	result$go_name = result$Term
	result$annotated = result$Annotated
	result$pval = test.result@score[result$go_id] # Rewrite p value to make numeric
	result$padj = p.adjust(result$pval,'fdr') # Calculate FDR
	
	# If the definition object exists has additional columns (besides pathway_id and pathway_name), tack them on
	# if (exists('definition',inherit=FALSE)) {
	if (!is.null(definition)) {
		definition.join = definition
		definition.join$go_id = definition.join$pathway_id
		definition.join$pathway_id = NULL
		definition.join$pathway_name = NULL
		result = merge(result,definition.join,by='go_id',all.x=TRUE)
	}
	
	result = result[order(result$pval),]
	rownames(result) = NULL
	
	required.columns = c('go_id','go_name','annotated','pval','padj')
	drop.columns = c('GO.ID','Term','Annotated','Significant','Expected')
	result[,c(required.columns,setdiff(colnames(result),c(required.columns,drop.columns)))]
}
