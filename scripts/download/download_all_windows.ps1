# Downloads public reference/GEO files into data\raw\<accession>\ and verifies SHA-256 where a reference hash is known.
# Safe: never overwrites a file whose hash differs (keeps the new copy as *.new). Logs to logs\download_windows.log
# Run by double-clicking RUN_DOWNLOADS.bat in the project folder (or: powershell -ExecutionPolicy Bypass -File scripts\download\download_all_windows.ps1)
$ErrorActionPreference = 'Continue'
$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
Set-Location $root
New-Item -ItemType Directory -Force -Path logs | Out-Null
$log = Join-Path $root 'logs\download_windows.log'
function Log($m){ $s = "$(Get-Date -Format s)  $m"; Write-Host $s; Add-Content -Path $log -Value $s }
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$ProgressPreference = 'SilentlyContinue'
$items = Import-Csv -Delimiter "`t" -Path (Join-Path $root 'data\download_list.tsv')
foreach ($it in $items) {
  $dir = Join-Path $root ("data\raw\" + $it.accession)
  New-Item -ItemType Directory -Force -Path $dir | Out-Null
  $out = Join-Path $dir $it.file_name
  if (Test-Path $out) {
    $h = (Get-FileHash -Algorithm SHA256 $out).Hash.ToLower()
    if ($it.sha256 -and $h -eq $it.sha256) { Log "OK present  $($it.file_name)"; continue }
    if (-not $it.sha256) { Log "present (no reference hash) $($it.file_name) sha256=$h"; continue }
    $tmp = "$out.new"
  } else { $tmp = "$out.part" }
  Log "GET $($it.url)"
  try { Invoke-WebRequest -Uri $it.url -OutFile $tmp -UseBasicParsing -TimeoutSec 3600 }
  catch { Log "  FAILED $($it.file_name): $($_.Exception.Message)"; continue }
  $h = (Get-FileHash -Algorithm SHA256 $tmp).Hash.ToLower()
  $size = (Get-Item $tmp).Length
  if ($it.sha256 -and $h -ne $it.sha256) { Log "  CHECKSUM MISMATCH $($it.file_name) got $h expected $($it.sha256); kept as $tmp"; continue }
  if ($tmp.EndsWith('.part')) { Move-Item -Force $tmp $out }
  Log "  verified $($it.file_name) bytes=$size sha256=$h"
}
Log "DONE"
