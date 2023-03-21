#' Balancing intercepts: update a definitions table to (approximately) control marginal means
#'
#' @param def Definition data.table to be modified
#' @param new_means Named numeric vector. Names are variables whose formula to modify; numbers are desired means.
#' @param omit Character vector. Names of variables to ignore (see `default`)
#' @param default Numeric. Desired marginal mean for all variables not mentioned in `new_means`. Default 0.
#' @param round_digits Integer, default 3. Number of digits to round coefficients (trust me, you want this - if you don't round you will get a gigantic number of digits and the table will be unreadible.)
#'
#' @return A data.table that is a data definition as in `simstudy`.
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
#' def = pt_def(formulas, dists, links=default_link(dists),
#'        variances, periods = 3, u_coef = 1)
#'
#' balance_intercepts(def,
#'                    new_means = c(a1_conditional =0.3, a2_conditional = 0.3, y0=1, y1=1.5, y2=-1),
#'                    omit=c('a1','a2'),
#'                    default = 0)
balance_intercepts = function(def, new_means, omit, default=0, round_digits=3) {

  get_balancing_intercept = function(formula, mean, var_means, link) {

    if(!grepl('\\+', formula) & link=='identity') {
      return(mean) #in this case the formula does not have covariates, just return the mean
    }
    if(!grepl('\\+', formula) & link=='logit') return(plogis(mean))

    mean_replacer = magrittr::set_names(as.character(var_means), names(var_means))

    formula_nointercept =  stringr::str_split(formula, '\\+', n=2, simplify = TRUE)[,2]
    formula_nointercept_with_means = stringr::str_replace_all(formula_nointercept, mean_replacer)

    balancer = if(link=='identity') mean else if(link=='logit') -log(1/mean - 1) else stop('link must be "logit" or "identity"')

    new_intercept = balancer - eval(parse(text = formula_nointercept_with_means))

    return( if(is.infinite(round_digits)) new_intercept else round(new_intercept, round_digits) )

  }

  replace_intercept = function(formula, new_intercept) {
    if(!grepl('\\+', formula)) {
      return(new_intercept)
    } else {
      split_formula = stringr::str_split(formula, pattern='\\+', n=2, simplify=TRUE)
      split_formula[,1] = new_intercept
      return(paste(split_formula, collapse='+') )
    }
  }

  indices_keep = if(!missing(omit)) seq_along(def$varname)[-match(omit, def$varname)] else seq_along(def$varname)
  varnames = def$varname[indices_keep]
  formulas = def$formula[indices_keep] %>% magrittr::set_names(varnames)
  links = def$link[indices_keep] %>% magrittr::set_names(varnames)

  #fill in defaults
  new_mean_vector = rep(default, length(def$varname)) %>% magrittr::set_names(def$varname)
  new_mean_vector[names(new_means)] = new_means

  #currently this will not get the means exactly right for anything with ak in it
  #it will just treat those as zeros. Probably not far off and good enough for our purposes.
  balancing_intercepts = sapply(varnames, function(v) get_balancing_intercept(formulas[v],
                                                                              mean = new_mean_vector[v],
                                                                              var_means = new_mean_vector,
                                                                              link = links[v])) %>%
    magrittr::set_names(varnames)



  new_formulas = sapply(varnames, function(v) replace_intercept(formulas[v], balancing_intercepts[v])) %>%
    magrittr::set_names(varnames)

  def$formula[indices_keep] = new_formulas
  def
}

#' Intervene on a data definition
#'
#' @param def Definition data.table to be modified
#' @param varnames Character vector. Names of variables whose value the intervention will modify.
#' @param values Numeric vector, same length as `varnames`. Values those variables will take under the intervention.
#'
#' @return A data.table that is an updated data definitions table.
#' @export
#'
#' @examples
#' def = pt_def(formulas = list(w1 = 'w1{t-1} + a{t-1}',
#'                              a  = 'w1{t} + u',
#'                              y  = 'w1{t} + a{t} + u'),
#'              dists = c('normal','binary','normal'),
#'              variances=c(1,1,1), periods = 3, u_coef = 1)
#'
#' def0 = intervene(def, varnames = c('a1','a2'), values = c(0,0))
#'
intervene = function(def, varnames, values) {
  stopifnot(length(varnames) == length(values))
  for(i in seq_along(varnames)) {
    def = updateDef(def, changevar = varnames[i], newformula = values[i], newdist = 'nonrandom')
  }
  return(def)
}
