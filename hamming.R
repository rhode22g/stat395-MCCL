# function copied from johanndejong.wordpress.com

hamming <- function(X) {
  D <- (1 - X) %*% t(X)
  D + t(D)
}