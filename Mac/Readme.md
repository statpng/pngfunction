## Sublime text 3
### [Sublime text setting](https://lazyren.github.io/devlog/sublime-text-setting.html)


## R
### [Install R]()
### [Solving the local warning problem](https://stackoverflow.com/questions/9689104/installing-r-on-mac-warning-messages-setting-lc-ctype-failed-using-c)
```
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C"
2: Setting LC_COLLATE failed, using "C"
3: Setting LC_TIME failed, using "C"
4: Setting LC_MESSAGES failed, using "C"
5: Setting LC_PAPER failed, using "C"
[R.app GUI 1.50 (6126) x86_64-apple-darwin9.8.0]

WARNING: You're using a non-UTF8 locale, therefore only ASCII characters will work. Please read R for Mac OS X FAQ (see Help) section 9 and adjust your system preferences accordingly.
```


### Running R on Terminal 
``` bash
export PATH="/Library/Frameworks/R.framework/Resources/bin:$PATH"
```
