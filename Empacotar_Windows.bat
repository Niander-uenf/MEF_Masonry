@echo off
setlocal
cd /d "%~dp0"
title MEF_Masonry - gerar executavel
echo Gerando o executavel do MEF_Masonry (dist\MEF_Masonry). Leva alguns minutos.
echo.
set "PY="
where py >nul 2>&1
if not errorlevel 1 (
  py --version
  if not errorlevel 1 set "PY=py"
)
if not defined PY (
  where python >nul 2>&1
  if not errorlevel 1 (
    python --version
    if not errorlevel 1 set "PY=python"
  )
)
if not defined PY (
  echo Nao encontrei o Python. Instale o Python 3.13 e rode este arquivo de novo.
  pause
  goto fim
)
echo Conferindo as bibliotecas...
%PY% -m pip install -r requirements.txt
%PY% empacotar.py
if errorlevel 1 (
  echo.
  echo A geracao falhou. Copie as mensagens acima e envie para o suporte.
) else (
  echo.
  echo Pronto. Copie a pasta dist\MEF_Masonry inteira para distribuir o programa.
)
pause
:fim
endlocal
