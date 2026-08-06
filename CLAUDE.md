# Project Configuration

## Kyber Implementation Directory
All Kyber work is in `elmo/projects/Examples/Kyber/`. Do NOT modify files outside of `/Users/sanjiv/Documents/Research/python-elmo/elmo/projects/Examples/Kyber/`.

## Coding Standards
- C90 compatible (ARM-none-eabi, Cortex-M0 target)
- No floating point (msoft-float)
- No dynamic memory allocation
- No standard library functions beyond stdint, string

## Build System
- Uses ARM toolchain: gcc, objcopy, objdump
- Add new files to SOURCES, HEADERS, OBJECTS in Makefile
- Build: make (outputs project.elf, project.bin, project.list)
