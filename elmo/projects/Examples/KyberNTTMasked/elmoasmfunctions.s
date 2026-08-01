##
## University of Bristol – Open Access Software Licence
## Copyright (c) 2016, The University of Bristol, a chartered
## corporation having Royal Charter number RC000648 and a charity
## (number X1121) and its place of administration being at Senate
## House, Tyndall Avenue, Bristol, BS8 1TH, United Kingdom.
## All rights reserved
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions
## are met:
##
## 1. Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
##
## 2. Redistributions in binary form must reproduce the above
## copyright notice, this list of conditions and the following
## disclaimer in the documentation and/or other materials provided
## with the distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
## FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
## COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
## HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
## STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
## OF THE POSSIBILITY OF SUCH DAMAGE.
##
## Any use of the software for scientific publications or commercial
## purposes should be reported to the University of Bristol
## (OSI-notifications@bristol.ac.uk and quote reference 2668). This is
## for impact and usage monitoring purposes only.
##
## Enquiries about further applications and development opportunities
## are welcome. Please contact elisabeth.oswald@bristol.ac.uk
##

.syntax unified
.text
.thumb

.func starttrigger
.thumb_func
.global starttrigger

starttrigger:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    movs r5, #0x04
    eors r4, r5
    movs r5, #1
    str r5, [r4, #0]
    pop {r0-r7}
bx lr

.endfunc

.func endtrigger
.thumb_func
.global endtrigger

endtrigger:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    movs r5, #0x04
    eors r4, r5
    movs r5, #0
    str r5, [r4, #0]
    pop {r0-r7}
bx lr

.endfunc

.func readbyte
.thumb_func
.global readbyte

readbyte:
    push {r0-r7}
    movs r4, #0xE1
    lsls r4, #24
    ldr r1, [r4, #0]
    strb r1, [r0, #0]
    pop {r0-r7}

bx lr

.endfunc

.func LoadN
.thumb_func
.global LoadN

LoadN:
    push {r0-r7}
    movs r4, #0xE1
    lsls r4, #24
    movs r5, #0x10
    eors r4, r5
    ldr r1, [r4, #0]
    str r1, [r0, #0]
    pop {r0-r7}

bx lr

.endfunc

.func randbyte
.thumb_func
.global randbyte

randbyte:
    push {r0-r7}
    movs r4, #0xE1
    lsls r4, #24
    movs r5, #0x04
    eors r4, r5
    ldr r1, [r4, #0]
    strb r1, [r0, #0]
    pop {r0-r7}

bx lr

.endfunc

.func getstart
.thumb_func
.global getstart

getstart:
    push {r0-r7}
    movs r4, #0xE1
    lsls r4, #24
    movs r5, #0x08
    eors r4, r5
    ldr r1, [r4, #0]
    str r1, [r0, #0]
    pop {r0-r7}

bx lr

.endfunc


.func getruncount
.thumb_func
.global getruncount

getruncount:
    push {r0-r7}
    movs r4, #0xE1
    lsls r4, #24
    movs r5, #0x0C
    eors r4, r5
    ldr r1, [r4, #0]
    str r1, [r0, #0]
    pop {r0-r7}

bx lr

.endfunc



.func printbyte
.thumb_func
.global printbyte

printbyte:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    ldrb r5, [r0]
    str r5, [r4]
    pop {r0-r7}
bx lr

.endfunc

.func endprogram
.thumb_func
.global endprogram

endprogram:
    push {r0-r7}
    movs r4, #0xF0
    lsls r4, #24
    movs r5, #0

    str r5, [r4]
    pop {r0-r7}
bx lr

.endfunc

.func initialisemaskflow
.thumb_func
.global initialisemaskflow

# Takes address of key as input (r0)
initialisemaskflow:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    movs r5, #0x40
    eors r4, r5
    str r0, [r4]
    pop {r0-r7}
bx lr

.endfunc

.func resetmaskflow
.thumb_func
.global resetmaskflow

resetmaskflow:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    movs r5, #0x42
    eors r4, r5
    movs r5, #0
    str r5, [r4]
    pop {r0-r7}
bx lr

.endfunc

.func setmaskflowstart
.thumb_func
.global setmaskflowstart

# Takes r0 as start number

setmaskflowstart:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    movs r5, #0x44
    eors r4, r5
    str r0, [r4]
    pop {r0-r7}
bx lr

.endfunc

.func resetdatafile
.thumb_func
.global resetdatafile

resetdatafile:
    push {r0-r7}
    movs r4, #0xE0
    lsls r4, #24
    movs r5, #0x46
    eors r4, r5
    movs r5, #0
    str r5, [r4]
    pop {r0-r7}
bx lr

.endfunc
