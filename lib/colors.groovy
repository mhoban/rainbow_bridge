//Regular text

class colors {
  static public String black = "\033[0;30m"
  static public String red = "\033[0;31m"
  static public String green = "\033[0;32m"
  static public String yellow = "\033[0;33m"
  static public String blue = "\033[0;34m"
  static public String magenta = "\033[0;35m"
  static public String cyan = "\033[0;36m"
  static public String white = "\033[0;37m"

  static public String black(str) {
    return "${this.black}${str}${this.reset}"
  }
  static public String red(str) {
    return "${this.red}${str}${this.reset}"
  }
  static public String green(str) {
    return "${this.green}${str}${this.reset}"
  }
  static public String yellow(str) {
    return "${this.yellow}${str}${this.reset}"
  }
  static public String blue(str) {
    return "${this.blue}${str}${this.reset}"
  }
  static public String magenta(str) {
    return "${this.magenta}${str}${this.reset}"
  }
  static public String cyan(str) {
    return "${this.cyan}${str}${this.reset}"
  }
  static public String white(str) {
    return "${this.white}${str}${this.reset}"
  }

  //Regular bold text
  static public String bblack = "\033[1;30m"
  static public String bred = "\033[1;31m"
  static public String bgreen = "\033[1;32m"
  static public String byellow = "\033[1;33m"
  static public String bblue = "\033[1;34m"
  static public String bmagenta = "\033[1;35m"
  static public String bcyan = "\033[1;36m"
  static public String bwhite = "\033[1;37m"

  static public String bblack(str) {
    return "${this.bblack}${str}${this.reset}"
  }
  static public String bred(str) {
    return "${this.bred}${str}${this.reset}"
  }
  static public String bgreen(str) {
    return "${this.bgreen}${str}${this.reset}"
  }
  static public String byellow(str) {
    return "${this.byellow}${str}${this.reset}"
  }
  static public String bblue(str) {
    return "${this.bblue}${str}${this.reset}"
  }
  static public String bmagenta(str) {
    return "${this.bmagenta}${str}${this.reset}"
  }
  static public String bcyan(str) {
    return "${this.bcyan}${str}${this.reset}"
  }
  static public String bwhite(str) {
    return "${this.bwhite}${str}${this.reset}"
  }


  //Regular underline text
  static public String UBLK = "\033[4;30m"
  static public String URED = "\033[4;31m"
  static public String UGRN = "\033[4;32m"
  static public String UYEL = "\033[4;33m"
  static public String UBLU = "\033[4;34m"
  static public String UMAG = "\033[4;35m"
  static public String UCYN = "\033[4;36m"
  static public String UWHT = "\033[4;37m"

  //Regular background
  static public String black_bg = "\033[40m"
  static public String red_bg = "\033[41m"
  static public String green_bg = "\033[42m"
  static public String yellow_bg = "\033[43m"
  static public String blue_bg = "\033[44m"
  static public String magenta_bg = "\033[45m"
  static public String cyan_bg = "\033[46m"
  static public String white_bg = "\033[47m"

  //High intensty background 
  static public String BLKHB = "\033[0;100m"
  static public String REDHB = "\033[0;101m"
  static public String GRNHB = "\033[0;102m"
  static public String YELHB = "\033[0;103m"
  static public String BLUHB = "\033[0;104m"
  static public String MAGHB = "\033[0;105m"
  static public String CYNHB = "\033[0;106m"
  static public String WHTHB = "\033[0;107m"

  //High intensty text
  static public String HBLK = "\033[0;90m"
  static public String HRED = "\033[0;91m"
  static public String HGRN = "\033[0;92m"
  static public String HYEL = "\033[0;93m"
  static public String HBLU = "\033[0;94m"
  static public String HMAG = "\033[0;95m"
  static public String HCYN = "\033[0;96m"
  static public String HWHT = "\033[0;97m"

  //Bold high intensity text
  static public String BHBLK = "\033[1;90m"
  static public String BHRED = "\033[1;91m"
  static public String BHGRN = "\033[1;92m"
  static public String BHYEL = "\033[1;93m"
  static public String BHBLU = "\033[1;94m"
  static public String BHMAG = "\033[1;95m"
  static public String BHCYN = "\033[1;96m"
  static public String BHWHT = "\033[1;97m"

  //Reset
  static public String reset =  "\033[0m"
  static public String CRESET = "\033[0m"
  static public String COLOR_RESET = "\033[0m"
}

