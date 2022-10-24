# non-parametric-calibration-development
Taking TJHeaton's code for the same and developing it to make it suitable to run with package.

To run this code, the carbondate package must first be installed. Because the package is currently private a few extra steps are needed so that (the authorised users only) can install it.

```
library('dev-tools')

## set your user name and email:
usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")

## create a personal access token for authentication:
usethis::create_github_token() # This takes you to the page to create the personal access token

## set personal access token:
credentials::set_github_pat("YourPAT") # The token you've just created
```

After this is done, the final step is to run `install_github("TJHeaton/carbondate")`.

Any of the scripts in this directory can then be run.
