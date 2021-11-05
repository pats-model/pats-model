# Parcel Ascent Tracing System (PATS)

[![dependency status](https://deps.rs/repo/github/Quba1/PATS/status.svg)](https://deps.rs/repo/github/Quba1/PATS)
[![license](https://img.shields.io/github/license/Quba1/PATS)](https://choosealicense.com/licenses/gpl-3.0/)

Parcel Ascent Tracing System (PATS) is the numerical model for convective parcel ascent simulation in three-dimension developed as a Master's project at the School of Environemntal Sciences at the University of East Anglia.

**The model is currently under development.** It will predict how convective air parcel behaves in conditions much more realistic then those assumed by so-called parcel theory.

For the details on how to configure and run the model check the [documentation](https://quba1.github.io/PATS/).

## Unique features of PATS

PATS builds upon lessons learned when developing world's most popular numerical weather prediction models, such as WRF, ECMWF IFS, ICON or GFS. Therefore it has several unique features to make it easy-to-use without sacrificing the model performance or flexibility.

### Modern programming language

PATS is written entirely in [Rust](https://www.rust-lang.org/), a modern programming language comprising great performance comparable to C and C++, type and memory safety and fearless concurrency. Using Rust allows for easy model development now and painless maintenance in the future.

### Detailed documentation

PATS documentation aims to be as insightful as that of ECMWF IFS or ICON, but also user-firendly. It covers mathematical and physical background of model, details of code design and the instructions of model usage. Therefore it is helpful both for finding out how to use and how to contribute to the model.

### Strongly-typed configuration file

PATS uses [`.yaml`](https://en.wikipedia.org/wiki/YAML) file for model configuration deserialized by [`serde`](https://serde.rs/). Therefore the configuration file is not only checked if it contains all required variables, but if they are of correct type as well. And if not, the model will provide an actually useful error message.

### First-class error handling

Use of Rust allows for an easy error handling even in concurrent applications, and segmentation faults are extremly rare. Even if something goes terribly wrong, the Rust will [panic](https://doc.rust-lang.org/book/ch09-01-unrecoverable-errors-with-panic.html) providing at least some error information. All other errors will be displayed with details in `stdout`, for example:

```text
[2021-09-28T12:25:05.475Z ERROR pats] Model failed with error: Error while reading config.yaml: Cannot deserialize config.yaml: domain.ni: invalid type: floating point `100.2`, expected u16 at line 6 column 7
```

By default PATS will display all errors, warnings and infos from log. However, if that is not sufficient you can increase the logging level using `PATS_LOG_LEVEL=debug` environmental variable to turn on the display of insightful debug messages.

### Non-blocking I/O

PATS concurrency design allows for the model to not stop entirely when writing the output and reading input. The model  will lose some performance temporalily, as IO operations are CPU-intensive. But it will use as much of remaining resources as possible, so no time is wasted waiting for IO operation to finish.

### Conditional compilation of model schemes

Most of numerical models compile all available features and schemes. That results in long compilation times, making debugging difficult, and performance overhead. PATS uses conditional compilation to compile only selected schemes. When the model is used operationally it has best possible performance, and when tested it can be quickly recompiled.
