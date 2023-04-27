#!/usr/bin/env Rscript

setClass(
	Class = 'pathwayAnnotation',
	representation(
		source='character',
		organism='character',
		annotation='list',
		definition='data.frame',
		ensembl_genes_version='character',
		pathway_annotation_version='character'
	)
)

setMethod(
    f = 'show',
    signature = 'pathwayAnnotation',
    definition = function(object) {
	    cat(paste0(
                                                                    '           Annotation name : ',object@source,'\n',
            if (!is.null(object@pathway_annotation_version)) paste0('Pathway annotation version : ',object@pathway_annotation_version,'\n') else '',
            if (!is.null(object@ensembl_genes_version))      paste0('     Ensembl genes version : ',object@ensembl_genes_version,'\n') else '',
            if (!is.null(object@organism))                   paste0('                  Organism : ',object@organism,'\n\n') else '\n',
                                                                    '           Number of genes : ',format(length(object@annotation),big.mark=',',scientific=FALSE),'\n',
                                                                    '        Number of pathways : ',format(length(unique(unlist(object@annotation))),big.mark=',',scientific=FALSE),'\n'
    ))
})

setMethod(
    f = 'subset',
    signature = 'pathwayAnnotation',
    definition = function(x,expr, ...) {
		safe_eval = function (expr, envir, enclos = parent.env(envir)) {
			expr = eval(call('bquote', expr, enclos))
			eval(expr, envir, enclos)
		}
		i = safe_eval(substitute(expr), x@definition, parent.frame(), ...)
		if (all(i)) x else x[,i]
    }
)

setMethod(
    f = 'dim',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	c(length(x@annotation),length(unique(unlist(x@annotation))))
    }
)

setMethod(
    f = 'nrow',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	length(x@annotation)
    }
)

setMethod(
    f = 'length',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	length(x@annotation)
    }
)

setMethod(
    f = 'ncol',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	length(unique(unlist(x@annotation)))
    }
)

setMethod(
    f = 'names',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	names(x@annotation)
    }
)

setMethod(
    f = 'rownames',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	names(x@annotation)
    }
)

setMethod(
    f = 'colnames',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	as.character(x@definition$pathway_id)
    }
)

setMethod(
    f = 'dimnames',
    signature = 'pathwayAnnotation',
    definition = function(x) {
    	list(
    		names(x@annotation),
    		as.character(x@definition$pathway_id)
    	)
    }
)

setMethod(
    f = '[',
    signature = 'pathwayAnnotation',
    definition = function(x,i,j) {
    	# x = annotation object
    	# i = genes
    	# j = pathways
    	annotation = x@annotation
    	definition = x@definition

    	if (!missing(i)) {
    		if (class(i) != 'character') {
    			if (length(i) > nrow(x)) warning('numerical expression has ',length(i),' elements: only the first used')
    		}
    		if (class(i) != 'logical') i = unique(i)
    		annotation = annotation[i]
    	}
    	if (!missing(j)) {
    		if (class(j) != 'character') {
    			j = as.character(definition$pathway_id[j])
    		}
    		annotation = lapply(annotation,function(x) intersect(x,j))
    	}

    	createPathwayAnnotationObject(
    		database = x@source,
    		annotation = annotation,
    		definition = definition,
    		organism = x@organism,
    		ensembl.version = x@ensembl_genes_version,
    		annotation.version = x@pathway_annotation_version
    	)
    }
)

#' Create a pathwayAnnotation object
#'
#' This function binds several data structures together into a single convenient object,
#' including an annotation object (a named list named by gene ID and containing a vector
#' of associated pathway IDs) and a definitions object (a data frame containing the
#' columns "pathway_id" and "pathway_name".
#'
#' @param database Name of the annotation source.
#' @param annotation Named list where each list item is named by gene ID and contains vector of pathway IDs.
#' @param definition Data frame containing the columns `pathway_id`, `pathway_name`, and optionally additional columns.
#' @param organism Name of the organism. It is advisable to use ENSEMBL codes, which can be retrieved using `fetchOrganisms('ENSEMBL')`.
#' @param ensembl.version ENSEMBL version
#' @param annotation.version Annotation version
#' @return An object of pathwayAnnotation class
#' @export
createPathwayAnnotationObject = function(
	database,
	annotation,
	definition,
	organism = NULL,
	ensembl.version = NULL,
	annotation.version = NULL
) {
	if (class(annotation) != 'list' || is.null(names(annotation))) stop('annotation must be a named list object.')
	if (!'data.frame' %in% class(definition)) stop('definition must be a data.frame object.')
	if (!all(c('pathway_id','pathway_name') %in% colnames(definition))) stop('definition must include column names "pathway_id" and "pathway_name".')
	if (length(database) > 1) warning('database name exceeds a length of 1. Retaining only the first item.')

	annotation = annotation[as.logical(unlist(lapply(annotation,length)))]

	if (!identical(names(definition)[1:2],c('pathway_id','pathway_name'))) {
		if (ncol(definition) == 2) {
			definition = definition[c('pathway_id','pathway_name')]
		} else {
			definition = definition[c('pathway_id','pathway_name',setdiff(colnames(definition),c('pathway_id','pathway_name')))]
		}
	}

	# If any pathways are in the annotation but not definition, add them to the definition data frame
	pathways.to.add = setdiff(definition$pathway_id,unlist(annotation))
	if (length(pathways.to.add)) {
		pathways.to.add = do.call(data.frame,c(list(pathways.to.add),lapply(rep(NA,ncol(definition)-1),function(x) x)))
		names(pathways.to.add) = names(definition)
		definition = rbind(definition,pathways.to.add)
	}

	definition = subset(definition,pathway_id %in% unlist(annotation))
	rownames(definition) = NULL

	new(
		'pathwayAnnotation',
		source = database[1],
		organism = organism[1],
		annotation = annotation,
		definition = definition,
		ensembl_genes_version = ensembl.version[1],
		pathway_annotation_version = annotation.version[1]
	)
}

#' Get definitions from an annotation object.
#'
#' This function extracts the definitions data frame from a pathwayAnnotation object.
#'
#' @param annotation An object of pathwayAnnotation class.
#' @return A data frame with pathway metadata
#' @export
getDefinitions = function(annotation) {
	annotation@definition
}

#' Get genes version from an annotation object.
#'
#' This function extracts the ENSEMBL genes version from a pathwayAnnotation object.
#'
#' @param annotation An object of pathwayAnnotation class.
#' @return A string with the ENSEMBL genes version
#' @export
getGenesVersion = function(annotation) {
	annotation@ensembl_genes_version
}

#' Get pathway annotations version from an annotation object.
#'
#' This function extracts the pathway annotations version from a pathwayAnnotation object.
#'
#' @param annotation An object of pathwayAnnotation class.
#' @return A string with the pathway annotations version
#' @export
getAnnotationsVersion = function(annotation) {
	annotation@pathway_annotation_version
}

#' Get annotations list from an annotation object.
#'
#' This function extracts the annotations list from a pathwayAnnotation object.
#'
#' @param annotation An object of pathwayAnnotation class.
#' @return A named list with gene-to-pathway associations
#' @export
getAnnotationsList = function(annotation) {
	annotation@annotation
}
