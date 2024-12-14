# How to push changes as a developer.

## Generating/Updating environment
### Windows ###
```powershell
conda env export | Select-String -NotMatch "^prefix: " | Out-File -Encoding utf8 environment_windows.yaml
(Get-Content environment_windows.yaml) -replace '=[^=]+', '' | Set-Content environment.yaml
pip freeze > .\docs\requirements.txt
```

### Linux ###
```bash
conda env export | grep -v "^prefix: " | sed -E 's/(=|==)[^=]+//g' > environment_linux.yml
pip freeze > docs\requirements.txt
```
