#' Create a definitions table for a `simstudy` simulation design under parallel trends
#'
#' @param formulas named list.
#' @param dists character vector.
#' @param links character vector.
#' @param variances  numeric vector.
#' @param periods integer.
#' @param u_coef numeric vector of length 1 or `length(formulas)`.
#'
#' @return A data.table that is a data definitions table for simstudy.
#' @export
#'
#' @examples
#' formulas = list(
#'   w1 = 'w1{t-1} + a{t-1}', #note no intercepts - these will be added automatically
#'   w2 = 'w1{t}   + w2{t-1} + a{t-1}',
#'   w3 = 'w1{t}   + w2{t}   + w3{t-1} + a{t-1}',
#'   a  = 'w1{t}   + cos(w2{t})   + I(w3{t}^2)  + u',
#'   y  = 'sin(w1{t})   + w2{t}*w3{t}   + a{t} + u'
#' )
#' dists =  c(
#'  w1 = 'normal',
#'  w2 = 'normal',
#'  w3 = 'normal',
#'  a  = 'binary',
#'  y  = 'normal'
#' )
#' variances = c(
#'   w1 = .1,
#'   w2 = .1,
#'   w3=  .1,
#'   a = 1, #this has no effect
#'   y = 0.1
#' )
#' pt_def(formulas, dists, links=default_link(dists),
#'        variances, periods = 3, u_coef = 1)
#'
#'
pt_def = function(formulas, dists, links = default_link(dists), variances = rep(1, length(formulas)), periods, u_coef = 1) {

  stopifnot(check_equal_length(formulas, dists, links, variances))

  varnames <- names(dists) <- names(links) <- names(variances) <- names(formulas)


  if(!check_varnames(varnames)) stop('variable names in formula argument must be exactly a, y, and w{i} where i are numbers')
  stopifnot(length(u_coef) %in% c(1, periods)) #allowing u_coef to vary for violations of parallel trends

  n_covs =  varnames %>% stringr::str_starts('w') %>% sum()

  #t=0 variables
  def = simstudy::defData(varname = 'u', dist = 'normal', formula = 0, variance = 1)
  for(j in 1:n_covs) {
    def = simstudy::defData(def,
                            varname = as.character(glue('w{j}{t}', .envir = list(t = 0, j=j))),
                            dist = dists[glue('w{j}')],
                            formula = glue_formula(formulas[[glue('w{j}')]], glue_vars = list(t=0)) %>% add_coefs_to_formula(),
                            variance = variances[glue('w{j}')])
  }

  def = simstudy::defData(def,
                varname = as.character(glue('a{t}', .envir = list(t = 0))),
                dist = 'nonrandom', formula = 0)
  def = simstudy::defData(def,
                varname = as.character(glue('y{t}', .envir = list(t=0))),
                dist = dists['y'],
                link = links['y'],
                formula = glue_formula(formulas$y, glue_vars = list(t=0)) %>% add_coefs_to_formula(u_coef = u_coef[1]),
                variance = variances[glue('y')])

  #t = 1:tau variables
  for(k in 1:(periods - 1)) {
    #w first
    for(varname in varnames)
    {
      def = simstudy::defData(def,
                    varname = paste0(varname, k, if(varname=='a') '_conditional'),#first generate what A{k} should be conditional on A{k-1}=0
                    formula = glue_formula(formulas[varname], glue_vars = list(t=k)) %>%
                      add_coefs_to_formula(u_coef = if(varname !='y') NULL else if(length(u_coef)==1) u_coef else u_coef[k + 1]),
                    dist = dists[varname],
                    link = links[varname],
                    variance = variances[varname])

      if(varname == 'a')  def=simstudy::defData(def,  #then A{k} marginal
                                      varname = paste0(varname, k),
                                      formula = glue_formula('0 + a{t-1} + (1-a{t-1})*a{t}_conditional', glue_vars = list(t=k)),
                                      dist = 'nonrandom')
    }
  }

  return(def)

}



#' Add (optionally randomly-generated) coefficients to a formula in a definition data.table.
#'
#' @param chr_formula Character vector, of the form `'x1+x2+x3'`
#' @param coefs Numeric vector, coefficients, length equal to number of terms (including intercept)
#' @param u_coef Numeric length 1, coefficient for the variable `u`, which we fix to be constant over time to ensure parallel trends
#'
#' @return Another character vector with coefficients added to formulas, of the form `'a + b*x1 + c*x2 + d*x3'`
#' @export
#'
#' @examples
#' add_coefs_to_formula('x1+x2+x3')
add_coefs_to_formula = function(chr_formula, coefs = NULL, u_coef = NULL){

  if(length(chr_formula) == 0) {

    return(generate_coefficients(p = 1))

  } else {

    terms = terms.glue_formula(chr_formula)

    if(is.null(coefs)) {
      coefs = generate_coefficients(p = length(terms) + 1) #+1 for an intercept
    }

    if(!is.null(u_coef)) {
      coefs[ which(terms == 'u') + 1] = u_coef
    }

    stopifnot(length(coefs) == length(terms) + 1)
    coefs_times_terms = paste0(coefs, c('', paste0('*', terms)))

    return(paste(coefs_times_terms, collapse=' + '))

  }

}

generate_coefficients = function(p) round(rnorm(p), 2)

#' Get default link for a distribution
#'
#' @param dist chr. Currently 'normal', 'binomial', or 'binary'.
#'
#' @return Character vector containing link names.
#' @export
#'
#' @examples
#' default_link(c('normal','binary'))
default_link = function(dist) {
  dplyr::case_when(dist == 'normal' ~ 'identity',
                                               dist %in% c('binomial','binary') ~ 'logit',
                                               TRUE ~ 'identity')
}

check_varnames = function(varnames) {
  a_varnames = varnames[stringr::str_starts(varnames, 'a')]
  y_varnames = varnames[stringr::str_starts(varnames, 'y')]
  w_varnames = varnames[stringr::str_starts(varnames, 'w')]

  #check that all variable names start with w, a, or y
  starts_with_way = varnames %>% stringr::str_starts('w|a|y') %>% all()

  #contains all w, a, and y
  contains_way = 0 < length(a_varnames)*length(y_varnames)*length(w_varnames)

  #check that w variables are of form w{n} where n is a number
  w_form_correct = w_varnames %>% stringr::str_match('[^w+[:digit:]]') %>% is.na() %>% all()

  #check that a and y are just a and y
  a_form_correct = a_varnames %>% stringr::str_match('[^a]') %>% is.na() %>% all()
  y_form_correct = y_varnames %>% stringr::str_match('[^y]') %>% is.na() %>% all()

  return (starts_with_way & contains_way & w_form_correct & a_form_correct & y_form_correct)
}
