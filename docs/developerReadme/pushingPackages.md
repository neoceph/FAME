# How to publish to pypi as a developer.

### Linux ###
**1.** Update the changelog on CHANGELOG.md.

**2.** Commit changes and create a git tag with appropriate version.

**3.** Remove old distribution in the directory `rm -rf dist/*`

**4.** Update the version on pyproject.toml and README.md PYPI version.

**5.** Build the distribution with `python -m build`

**6.** Make sure `twine` is installed. If not, install using `pip install twine`

**7.** Before uploading make sure you have `.pypirc` file in your home directory with the following content
```bash
[testpypi]
  username = __token__
  password = token-from-test-pypi-account
```
or 
```bash
[pypi]
  username = __token__
  password = token-from-pypi-account
```
**8.** Upload using twine to testpypi
```bash
twine upload --repository testpypi dist/* --verbose
```
or to pypi
```bash
twine upload --repository pypi dist/* --verbose
```
