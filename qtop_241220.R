rm(list=ls())
gc()

library(dplyr)
library(stringr)
library(tidyselect)
library(purrr)
library(shinythemes)

ui <- fluidPage(
  
  tags$h2(strong("NONMEM file processor : qTMDD to pTMDD"), style = "font-size: 30px;"),
  
  tags$style(HTML("
    body {
      background-color: #F6F6F6;}

      .progress-bar {
        background-color:#A3B3B9; /* 업로드 완료 바의 색상 */
      }
  ")),
  sidebarLayout(
    sidebarPanel(
      radioButtons("processoption", "Choose processing option:",
                   choices = list("TMDD approximation for receptor and FcRn" = "yes", "TMDD approximation only for the receptor" = "no"), selected = "no"),
      radioButtons("origindv", "Origin of DV:",
                   choices = list("Free drug concentration" = "cfree", "Total drug concentration" = "ctot"), selected = "ctot"),
      numericInput("cmt_tmdd", "In which compartment the drug-receptor TMDD happens : A(n) ", 3, min = 1, max = 100),
      conditionalPanel(
        condition = "input.processoption == 'yes'",
        numericInput("cmt_fcrn", "In which compartment the drug-FcRn TMDD happens : A(n)", 2, min = 1, max = 100)
      ),
      fileInput("modFile", "Upload a .mod File", accept = c(".mod")),
      actionButton("process_btn", "Process File"),
      
      downloadButton("downloadData", "Download Processed File")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions",
                 img(src = "instruction3.png", height = "616px", width = "1011px")
        ),
        tabPanel("Processed File",
                 conditionalPanel(
                   condition = "input.process_btn == 0",
                   h5("Upload the file, choose the options, and then click the process button!")
                 ),
                 conditionalPanel(
                   condition = "input.process_btn > 0",
                   h4("<Processed File Output>"),
                   verbatimTextOutput("output_content")
                 )
        )
      )
    )
  )
)

server <- function(input, output) {
  processedData <- reactiveVal()  # 반응형 값으로 초기화
  
  observeEvent(input$process_btn, {
    req(input$modFile)
    # 파일 읽기 및 주석 제거
    myfile <- readLines(input$modFile$datapath) %>% gsub(";.*", "", .) %>% trimws() %>% .[. != ""]
    # $DES 및 $ERROR 섹션 찾기
    des <- which(grepl("\\$DES", myfile))
    error <- which(grepl("\\$ERROR", myfile))
    theta <- which(grepl("\\$THETA", myfile))
    
    # Extract the lines between $DES and $ERROR (inclusive)
    before <- myfile[1:(des-1)]
    part_des <- myfile[(des+1):(error-1)] 
    part_error <- myfile[(error+1):(theta-1)]
    after <- myfile[theta:length(myfile)]
    
    part_dadt <- grep("^DADT", part_des, value = TRUE)
    part_pre <- part_des[!part_des %in% part_dadt]
    
    # ------ From the DAA line, define ctot, kss, rtot ... and then Cfree
    # fcrn recycling part
    DF <- part_pre[grepl("^\\s*([^=]+)\\s*=\\s*([^\\-]+)\\s*\\-\\s*([^\\-]+)\\s*\\-\\s*([^\\-]+)\\s*$", part_des)] 
    DF <- DF1[!is.na(DF1)]
    
    if (input$processoption == "yes") {
      DF1 <- grep("fcrn", DF, ignore.case = TRUE, value = TRUE)
      if (length(DF1) == 0) {
        # fcrn을 찾지 못한 경우 "error" 출력
        processedData("error")
        return()
      }
    } 
    
    pattern1 = "\\s*=\\s*|\\s*\\-\\s*" # 뭐 = 뭐 - 뭐 - 뭐, strsplit에 사용'
    pattern2 <- "\\s*.+\\s*\\*\\s*.+\\s*/\\s*\\(.+\\s*\\+\\s*.+\\)"    # 뭐*뭐/(뭐+뭐) 형태 찾기, 공백 포함
    # For the Ccpx,
    # 미리 정의된 경우,
    if(input$processoption == "no"){
      
      DF2 <- DF
      DAA2_def <- (strsplit(DF2, pattern1))[[1]]
      DAA2 <- DAA2_def[1]
      Ctot2 <- DAA2_def[2]
      DAA2_left <-  DAA2_def[3:4]
      DAA2_K <- grep("^k", DAA2_def, value=TRUE, ignore.case=TRUE)
      Rtot2 <- DAA2_left[!DAA2_left %in% DAA2_K]
      
      index <- which(part_pre == DF2)
      Cfree2 <- part_pre[index+1]
      cf2 <- sub("^(.*?)\\s*=.*$", "\\1", Cfree2)
      
      pattern2 <- "\\s*.+\\s*\\*\\s*.+\\s*/\\s*\\(.+\\s*\\+\\s*.+\\)"    # 뭐*뭐/(뭐+뭐) 형태 찾기, 공백 포함
      matching_items <- grep(pattern2, part_pre, value = TRUE)
      
      if(length(matching_items)== 1){
        
        index2 <- which(part_pre == DF2)
        cpx2 <- part_pre[index2+2]
        cpx2 <- sub("^(.*?)\\s*=.*$", "\\1", cpx2)
        
        ccpx2 <- paste0(cpx2, " = ", Rtot2, " * ", Ctot2, " / (", DAA2_K, " + ", Ctot2, " + ", Rtot2, ")")
        Af2 <- paste0(cf2, " = ", Ctot2, "-", cpx2)
        
        part_pre_decluded <- part_pre[-c(index2, index2+1)]
        pattern <- "\\b\\S+\\s*\\*\\s*\\S+\\s*/\\s*\\(\\s*(?:\\S+\\s*\\*\\s*)?\\S+\\s*\\+\\s*\\S+\\s*\\)"
        part_pre <- part_pre_decluded[grep(pattern, part_pre_decluded, invert=TRUE)]
        
        part_pre <- c(part_pre, ccpx2, Af2)
        
      }
      
      else {
        
        index2 <- which(part_pre == DF2)
        ccpx2 <- paste0("Ccpx2", " = ", Rtot2, " * ", Ctot2, " / (", DAA2_K, " + ", Ctot2, " + ", Rtot2, ")")
        Af2 <- paste0(cf2, " = ", Ctot2, "-", "Ccpx2")
        
        cmtpattern1 <- paste0("DADT(",input$cmt_tmdd,")")
        
        line_rtot <- grep(cmtpattern1, part_des, value = TRUE, fixed=TRUE)
        
        # for the part + Fcrn recycled
        modified_line <- str_replace_all(line_rtot, "A\\((\\d+)\\)", "A\\1")
        modified_line <- str_replace_all(modified_line, "\\(([^()]*?)\\+([^()]*?)\\)", "\\(\\1 plus \\2\\)")
        modified_line <- str_replace_all(modified_line, "\\(([^()]*?)\\-([^()]*?)\\)", "\\(\\1 minus \\2\\)")
        
        # =, -, +을 기준으로 그 다음 항목 추출
        pattern <- "(=|\\+|\\-)\\s*[^=\\+\\-]+"
        result <- unlist(str_extract_all(modified_line, pattern))
        result <- trimws(result)
        # 정규표현식 패턴 정의: -로 시작하고 "/(뭐 plus 뭐)" 구조를 포함하는 항목
        pattern2 <- "^-\\s*[^/]+/\\s*\\([^\\)]+\\s+plus\\s+[^\\)]+\\)"
        target <- grep(pattern2, result, value=TRUE)
        
        # 추출된 항목을 *와 /, 그리고 -를 기준으로 분리
        split_target <- unlist(strsplit(target, "[*/-]")) %>% trimws()
        pattern3 <- "\\(\\s*[A-Za-z0-9]+\\s*plus\\s*[A-Za-z0-9]+\\s*\\)"
        denom <- grep(pattern3, split_target, value=TRUE) 
        
        group <- c(Rtot2, cf2, denom)
        
        target_left <- split_target[!split_target %in% group] %>% trimws()
        target_left <- target_left[target_left != ""]
        
        term_cp <- paste(paste(target_left, collapse = "*"), "*", "Ccpx2",  sep = "")
        target_left <- result[!result %in% target]
        result <- paste0("DADT(",input$cmt_tmdd,")", paste(target_left, collapse = ""), "-", term_cp, sep = " ") 
        
        result <- gsub("minus", "-", result)
        result <- gsub("plus", "+", result)
        
        # 결과 출력
        part_dadt[input$cmt_tmdd] <- result
        
        part_pre_decluded <- part_pre[-c(index2, index2+1)]
        pattern <- "\\b\\S+\\s*\\*\\s*\\S+\\s*/\\s*\\(\\s*(?:\\S+\\s*\\*\\s*)?\\S+\\s*\\+\\s*\\S+\\s*\\)"
        part_pre <- part_pre_decluded[grep(pattern, part_pre_decluded, invert=TRUE)]
        
        part_pre <- c(part_pre, ccpx2, Af2)
      }
      
    } else {
      
      DF1 <- grep("fcrn", DF, ignore.case = TRUE, value = TRUE)
      DAA1_def <- (strsplit(DF1, pattern1))[[1]]
      DAA1 <- DAA1_def[1]
      Ctot1 <- DAA1_def[2]
      DAA1_left <-  DAA1_def[3:4]
      DAA1_K <- grep("^k", DAA1_def, value=TRUE, ignore.case=TRUE)
      Rtot1 <- DAA1_left[!DAA1_left %in% DAA1_K]
      
      index <- which(part_pre == DF1)
      free1 <- part_pre[index+1]
      cf1 <- sub("^(.*?)\\s*=.*$", "\\1", free1)
      
      # drug receptor part
      DF2 <- DF[!DF %in% DF1][1]
      DAA2_def <- (strsplit(DF2, pattern1))[[1]]
      DAA2 <- DAA2_def[1]
      Ctot2 <- DAA2_def[2]
      DAA2_left <-  DAA2_def[3:4]
      DAA2_K <- grep("^k", DAA2_def, value=TRUE, ignore.case=TRUE)
      Rtot2 <- DAA2_left[!DAA2_left %in% DAA2_K]
      
      index <- which(part_pre == DF2)
      Cfree2 <- part_pre[index+1]
      cf2 <- sub("^(.*?)\\s*=.*$", "\\1", Cfree2)
      
      index1 <- which(part_pre == DF1)
      index2 <- which(part_pre == DF2)
      
      
      matching_items <- grep(pattern2, part_pre, value = TRUE)
      
      if(length(matching_items)== 2){
        
        cpx1 <- part_pre[index1+2]
        cpx1 <- sub("^(.*?)\\s*=.*$", "\\1", cpx1)
        cpx2 <- part_pre[index2+2]
        cpx2 <- sub("^(.*?)\\s*=.*$", "\\1", cpx2)
        
        ccpx1 <- paste0(cpx1, " = ", Rtot1, " * ", Ctot1, " / (", DAA1_K, " + ", Ctot1, " + ", Rtot1, ")")
        Af1 <- paste0(cf1, " = ", Ctot1, "-", cpx1)
        ccpx2 <- paste0(cpx2, " = ", Rtot2, " * ", Ctot2, " / (", DAA2_K, " + ", Ctot2, " + ", Rtot2, ")")
        Af2 <- paste0(cf2, " = ", Ctot2, "-", cpx2)
        
        part_pre[index+1]
        part_pre_decluded <- part_pre[-c(index1, index1 + 1, index2, index2+1)]
        pattern <- "\\b\\S+\\s*\\*\\s*\\S+\\s*/\\s*\\(\\s*(?:\\S+\\s*\\*\\s*)?\\S+\\s*\\+\\s*\\S+\\s*\\)"
        part_pre <- part_pre_decluded[grep(pattern, part_pre_decluded, invert=TRUE)]
        
        part_pre <- c(part_pre, ccpx1, Af1, ccpx2, Af2)
        
      } else {
        # for the predefined part,
        
        ccpx1 <- paste0("Ccpx1", " = ", Rtot1, " * ", Ctot1, " / (", DAA1_K, " + ", Ctot1, " + ", Rtot1, ")")
        Af1 <- paste0(cf1, " = ", Ctot1, "-", "Ccpx1")
        ccpx2 <- paste0("Ccpx2", " = ", Rtot2, " * ", Ctot2, " / (", DAA2_K, " + ", Ctot2, " + ", Rtot2, ")")
        Af2 <- paste0(cf2, " = ", Ctot2, "-", "Ccpx2")
        
        cmtpattern1 <- paste0("DADT(",input$cmt_tmdd,")")
        cmtpattern2 <- paste0("DADT(",input$cmt_fcrn,")")
        
        line_rtot <- grep(cmtpattern1, part_des, value = TRUE, fixed=TRUE)
        line_fcrn <- grep(cmtpattern2, part_des, value = TRUE, fixed=TRUE)
        
        # 정규 표현식을 사용하여 괄호로 묶인 항목과 일반 항목을 분리
        # 괄호 포함된 항목은 하나의 단위로 추출하고, 괄호 없는 항목은 `=`, `+`, `-`로 분리 modified_line <- str_replace_all(line, "A\\((\\d+)\\)", "A\\1")
        modified_line <- str_replace_all(line_fcrn, "A\\((\\d+)\\)", "A\\1")
        modified_line <- str_replace_all(modified_line, "\\(([^()]*?)\\+([^()]*?)\\)", "\\(\\1 plus \\2\\)")
        modified_line <- str_replace_all(modified_line, "\\(([^()]*?)\\-([^()]*?)\\)", "\\(\\1 minus \\2\\)")
        
        # =, -, +을 기준으로 그 다음 항목 추출
        pattern <- "(=|\\+|\\-)\\s*[^=\\+\\-]+"
        result <- unlist(str_extract_all(modified_line, pattern))
        result <- trimws(result)
        # 정규표현식 패턴 정의: -로 시작하고 "/(뭐 plus 뭐)" 구조를 포함하는 항목
        pattern2 <- "^-\\s*[^/]+/\\s*\\([^\\)]+\\s+plus\\s+[^\\)]+\\)"
        target <- grep(pattern2, result, value=TRUE)
        
        # 추출된 항목을 *와 /, 그리고 -를 기준으로 분리
        split_target <- unlist(strsplit(target, "[*/-]")) %>% trimws()
        pattern3 <- "\\(\\s*[A-Za-z0-9]+\\s*plus\\s*[A-Za-z0-9]+\\s*\\)"
        denom <- grep(pattern3, split_target, value=TRUE) 
        
        group <- c(Rtot1, cf1, denom)
        
        target_left <- split_target[!split_target %in% group] %>% trimws()
        target_left <- target_left[target_left != ""]
        
        term_cp <- paste(paste(target_left, collapse = "*"), "*", "Ccpx1",  sep = "")
        target_left <- result[!result %in% target]
        result <- paste0("DADT(",input$cmt_fcrn,")", paste(target_left, collapse = ""), "-", term_cp, sep = " ") 
        
        result <- gsub("minus", "-", result)
        result <- gsub("plus", "+", result)
        
        # 결과 출력
        part_dadt[input$cmt_fcrn] <- result
        
        # for the part + Fcrn recycled
        modified_line <- str_replace_all(line_rtot, "A\\((\\d+)\\)", "A\\1")
        modified_line <- str_replace_all(modified_line, "\\(([^()]*?)\\+([^()]*?)\\)", "\\(\\1 plus \\2\\)")
        modified_line <- str_replace_all(modified_line, "\\(([^()]*?)\\-([^()]*?)\\)", "\\(\\1 minus \\2\\)")
        pattern <- "(=|\\+|\\-)\\s*[^=\\+\\-]+"
        result <- unlist(str_extract_all(modified_line, pattern))
        result <- trimws(result)
        
        pattern2.0 <- "^\\+\\s*[^/]+/\\s*\\([^\\)]+\\s+plus\\s+[^\\)]+\\)"
        target2.0 <- grep(pattern2.0, result, value=TRUE)
        
        target_left <- result[!result %in% target2.0] 
        
        result <- paste0("DADT(",input$cmt_fcrn,")", paste(target_left, collapse = ""), "+", term_cp, sep = " ") 
        
        # for the part - TMDD with receptor
        # 정규표현식 패턴 정의: -로 시작하고 "/(뭐 plus 뭐)" 구조를 포함하는 항목
        result <- unlist(str_extract_all(result, pattern))
        target <- grep(pattern2, result, value=TRUE)
        
        # 추출된 항목을 *와 /, 그리고 -를 기준으로 분리
        split_target <- unlist(strsplit(target, "[*/-]")) %>% trimws()
        pattern3 <- "\\(\\s*[A-Za-z0-9]+\\s*plus\\s*[A-Za-z0-9]+\\s*\\)"
        denom <- grep(pattern3, split_target, value=TRUE) 
        
        group <- c(Rtot2, cf2, denom)
        
        target_left <- split_target[!split_target %in% group] %>% trimws()
        target_left <- target_left[target_left != ""]
        
        term_cp <- paste(paste(target_left, collapse = "*"), "*", "Ccpx2",  sep = "")
        target_left <- result[!result %in% target]
        result <- paste0("DADT(",input$cmt_tmdd,")", paste(target_left, collapse = ""), "-", term_cp, sep = " ") 
        
        result <- gsub("minus", "-", result)
        result <- gsub("plus", "+", result)
        part_dadt[input$cmt_tmdd] <- result
        
        
        part_pre[index+1]
        part_pre_decluded <- part_pre[-c(index1, index1 + 1, index2, index2+1)]
        pattern <- "\\b\\S+\\s*\\*\\s*\\S+\\s*/\\s*\\(\\s*(?:\\S+\\s*\\*\\s*)?\\S+\\s*\\+\\s*\\S+\\s*\\)"
        part_pre <- part_pre_decluded[grep(pattern, part_pre_decluded, invert=TRUE)]
        
        part_pre <- c(part_pre, ccpx1, Af1, ccpx2, Af2)
      }
    }
    # part_error
    if (input$origindv == 'cfree'){
      # 뭐 = 뭐 - 뭐 - 뭐 구조 추출 후 분리
      pattern <- "\\b\\S+\\s*=\\s*\\S+\\s*-\\s*\\S+\\s*-\\s*\\S+"
      DAA_def <- part_error[grep(pattern, part_error)]
      index <- which(part_error == DAA_def)
      DAA_def <- (strsplit(DAA_def, pattern1))[[1]]
      
      DAA <- DAA_def[1]
      Ctot <- DAA_def[2]
      DAA_left <-  DAA_def[3:4]
      DAA_k <- grep("^k", DAA_def, value=TRUE, ignore.case=TRUE)
      Rtot <- DAA_left[!DAA_left %in% DAA_k]
      
      cf <- part_error[index+1]
      cf <- sub("^(.*?)\\s*=.*$", "\\1", cf)
      
      ccpx <- paste0("Ccpx", " = ", Rtot, " * ", Ctot, " / (", DAA_k, " + ", Ctot, " + ", Rtot, ")")
      Af <- paste0(cf, " = ", Ctot, "-", "Ccpx")
      
      left <- part_error[index-1]
      left2 <- part_error[(index+2):length(part_error)]
      
      part_error <- c(left, ccpx, Af, left2)}  else {part_error <- part_error}
    
    all <- c(before,"\n","$DES",part_pre,"\n",part_dadt,"\n","$ERROR",part_error,"\n",after)
    
    output_content <- paste(all, collapse = "\n")
    processedData(output_content)
  })
  
  # 처리된 내용을 텍스트로 표시
  output$output_content <- renderText({
    req(processedData())
    processedData()
  })
  
  # 파일 다운로드
  output$downloadData <-  downloadHandler(
    filename = function() {
      original_name <- tools::file_path_sans_ext(input$modFile$name)  # 업로드된 파일의 원래 이름 가져오기
      paste("processed_", original_name, ".mod", sep = "")
    },
    content = function(file) {
      writeLines(processedData(), file)
    }
  )
}


# 앱 실행
shinyApp(ui = ui, server = server)