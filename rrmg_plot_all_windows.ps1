[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$WorkDir,

    [string]$CodeDir = $PSScriptRoot,

    [string]$Rscript,

    [string]$OutputDir
)

$ErrorActionPreference = 'Stop'
$WorkDir = (Resolve-Path -LiteralPath $WorkDir).Path
$CodeDir = (Resolve-Path -LiteralPath $CodeDir).Path

if (-not $Rscript) {
    $command = Get-Command Rscript.exe -ErrorAction SilentlyContinue
    if ($command) {
        $Rscript = $command.Source
    } else {
        $candidates = @(Get-ChildItem -Path 'C:\Program Files\R\R-*\bin\Rscript.exe' -File -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending)
        if ($candidates.Count -eq 1) {
            $Rscript = $candidates[0].FullName
        } elseif ($candidates.Count -gt 1) {
            throw "Multiple Rscript.exe installations found; provide -Rscript explicitly."
        } else {
            throw "Rscript.exe not found on PATH or under C:\Program Files\R."
        }
    }
}
if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) {
    throw "Rscript.exe not found: $Rscript"
}

$complete = Join-Path $WorkDir 'figure_data\COMPLETE.json'
$alignmentComplete = Join-Path $WorkDir 'figure_data\CYP51A_alignment\COMPLETE.json'
foreach ($path in @($complete, $alignmentComplete)) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Incomplete figure-data bundle: $path"
    }
}

if (-not $OutputDir) {
    $OutputDir = Join-Path $WorkDir 'figure_output'
}
$previewOut = Join-Path $OutputDir 'local_preview'
$alignmentOut = Join-Path $OutputDir 'alignment_sites'
$previewScript = Join-Path $CodeDir 'rrmg_plot_preview.R'
$alignmentScript = Join-Path $CodeDir 'rrmg_plot_alignment.R'

foreach ($path in @($previewScript, $alignmentScript)) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "R plotting script not found: $path"
    }
}

Write-Host "[R] Overview figures -> $previewOut"
& $Rscript $previewScript $WorkDir $previewOut
if ($LASTEXITCODE -ne 0) { throw "rrmg_plot_preview.R failed with exit code $LASTEXITCODE" }

Write-Host "[R] Alignment figures -> $alignmentOut"
& $Rscript $alignmentScript (Join-Path $WorkDir 'figure_data\CYP51A_alignment') $alignmentOut
if ($LASTEXITCODE -ne 0) { throw "rrmg_plot_alignment.R failed with exit code $LASTEXITCODE" }

Write-Host "[DONE] RRMG Windows figures: $OutputDir"
