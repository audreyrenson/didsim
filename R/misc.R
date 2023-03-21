check_equal_length = function(...) {
  list_vecs = list(...)
  length(unique(sapply(list_vecs, length))) == 1
}


glue_formula = function(string, glue_vars) {

  trimmed_string = stringr::str_remove_all(string, '[:space:]')
  glued_trimmed_string = glue(trimmed_string, .envir = glue_vars)

  remove_negatives(glued_trimmed_string)
}

remove_negatives = function(glued_string) {
  #this is to get rid of references to variables before time 0
  result = glued_string %>%
    stringr::str_remove_all('\\+[:alpha:]+-[:digit:]+') %>%#first remove those starting with +
    stringr::str_remove_all('\\*[:alpha:]+-[:digit:]+') %>%#or a *
    stringr::str_remove_all('\\:[:alpha:]+-[:digit:]+') %>%#etc.
    stringr::str_remove_all('[:alpha:]+-[:digit:]+') %>%
    stringr::str_remove_all('\\+[:alpha:]+[:digit:]+-[:digit:]+') %>%#first remove those starting with +
    stringr::str_remove_all('\\*[:alpha:]+[:digit:]+-[:digit:]+') %>%#or a *
    stringr::str_remove_all('\\:[:alpha:]+[:digit:]+-[:digit:]+') %>%#etc.
    stringr::str_remove_all('[:alpha:]+[:digit:]+-[:digit:]+')

  if(result == '') return(character(0)) else return(result)
}

terms.glue_formula = function(glued_string) {
  glued_string %>%
    paste0('~', .) %>%
    as.formula() %>%
    terms() %>%
    attr('term.labels') %>%
    stringr::str_replace_all(':', '*')
}
