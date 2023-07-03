# Load required libraries
library(shiny)
library(ssh)
library(tidytable)
library(ggplot2)
library(googlesheets4)
library(stringr)
library(forcats)



# Define the UI
ui <- fluidPage(
  # Main app page
  mainPanel(

    # Plots
    
    fluidRow(
      column(width=8, class = "well", 
             plotOutput("plot1",height = "300px"),
             p(em("Note:"),"MpileupDone is the final status of a run before AHMM")
      
    )),br(),
    fluidRow(
      column(width=8, class = "well", 
             plotOutput("plot2",height = "300px"),
             p(em("Note:"),"VCF should be 1 file per individual, not sequencing file, while the initial mapping is one file per sample per run")
      
    )),br(),
    fluidRow(
      column(width=8, class = "well", 
             plotOutput("plot3",height = "300px"),
             p(em("Note:"),"here is a quick check if the server is still busy")
      
    ))
  )
)

get_QC_df<-function(ssh_session,fx_directory){
  out <- ssh_exec_internal(ssh_session, paste0("cat ",fx_directory,"QC.tsv"))
  return(fread(text=rawToChar(out$stdout)))
}
get_status_df<-function(ssh_session,fx_directory){
  out <- ssh_exec_internal(ssh_session, paste0("cat ",fx_directory,"status.tsv"))
  return(fread(text=rawToChar(out$stdout)))
}
get_ahmm_files<-function(ssh_session,fx_directory){
  ahmm_files <- ssh_exec_internal(ssh_session, paste0("ls ",fx_directory,"*/*vcf.gz"))
  ahmm_files <- rawToChar(ahmm_files$stdout)
  ahmm_files <- strsplit(ahmm_files, "\n")[[1]]
  
  ahmm_files_mtime=c()
  # break it down into few commands
  ahmm_seq<-round(seq(1, length(ahmm_files), length.out=round(length(ahmm_files)/100)))
  if(length(ahmm_seq)<2){
    ahmm_seq<-c(1, length(ahmm_files))
  }
  print("Checking VCF file timestamps...")
  ahmm_seq[length(ahmm_seq)]<-ahmm_seq[length(ahmm_seq)]+1 # solution for the way I iterate
  for(i in 2:length(ahmm_seq)){
    af<-ahmm_files[ahmm_seq[i-1]:(ahmm_seq[i]-1)]
    af<-paste(af, collapse=" ")
    tmp_mtime<-ssh_exec_internal(ssh_session, paste0("stat -c '%Y' ",af))
    tmp_mtime<-rawToChar(tmp_mtime$stdout)
    tmp_mtime<-strsplit(tmp_mtime, "\n")[[1]]
    tmp_mtime<-as.numeric(tmp_mtime)
    ahmm_files_mtime<-c(ahmm_files_mtime, tmp_mtime)
    print(paste0(round(100*(ahmm_seq[i]/ahmm_seq[length(ahmm_seq)]),2),"%"))
  }
  print(length(ahmm_files))
  print(length(ahmm_files_mtime))
  return(data.table(file=ahmm_files, mtime=ahmm_files_mtime))
}
get_slurm_jobs<-function(ssh_session,username){
  today_slurm<-ssh_exec_internal(ssh_session, paste0("sacct -u ",username))
  today_slurm<-rawToChar(today_slurm$stdout)
  today_slurm<-strsplit(today_slurm, "\n")[[1]]
  today_slurm<-today_slurm[-1*c(1,2)]
  today_slurm_df<-data.table()
  for(ts in today_slurm){
    ts_split<-strsplit(ts," ")[[1]]
    ts_split<-ts_split[ts_split!=""]
    tmp_df<-data.table(jobname=ts_split[2],
                       status=ts_split[6])
    today_slurm_df<-bind_rows(today_slurm_df, tmp_df)
  }
  
  Fx_jobs<-c("BCL2Fastq",
             "DemultiFastq",
             "fx_trim",
             "map-fx",
             "QC-bam-fx",
             "mpileup-fx",
             "fx-sexing-prep",
             "mpileup-cors",
             "AHMM")
  Fx_jobs_plot<-c("BCL2Fastq",
                  "Demultiplex",
                  "Trimming",
                  "Mapping",
                  "Get QC data",
                  "Mpileup",
                  "Prep sexing data",
                  "Correlate sample\nallele counts\nfor merging",
                  "Run AHMM &\nconvert to VCF")
  
  today_slurm_df<-filter(today_slurm_df, jobname %in% Fx_jobs)
  recode_vector <- setNames(Fx_jobs, Fx_jobs_plot)
  
  today_slurm_df <- today_slurm_df %>%
    mutate(jobname = fct_recode(jobname, !!!recode_vector))%>%
    mutate(jobname = fct_relevel(jobname, Fx_jobs_plot)) # rename levels
  
  return(today_slurm_df)
}

# Define the server function
server <- function(input, output, session) {
  # Placeholder for username and directory
  fx_directory <- reactiveVal(NULL)
  username <- reactiveVal(NULL)
  jobs_username <- reactiveVal(NULL)
  password <- reactiveVal(NULL)
  ssh_session <- reactiveVal(NULL)
  status_df <- reactiveVal(NULL)
  QC_df <- reactiveVal(NULL)
  ahmm_files <- reactiveVal(NULL)
  files_df_fx<- reactiveVal(NULL)
  today_slurm_df<-reactiveVal(NULL)
  
  # Handle login
  observeEvent(input$login, {
    req(input$jobs_username, input$username, input$password, input$fx_directory)  # Ensure inputs are provided
    fx_directory(input$fx_directory)
    jobs_username(input$jobs_username)
    username(input$username)
    password(input$password)
    ssh_session(ssh_connect(paste0(username(),"@server"), passwd = password()))
    gs4_auth(email = paste0(username(),"@server"),cache="~/gcache")
    files_df_fx(read_sheet("sheet_id", sheet = 3))
    removeModal(session)
  })
  
  # Define placeholder plots
  output$plot1 <- renderPlot({
    req(ssh_session()) 
    QC_df(get_QC_df(ssh_session(),fx_directory()))
    status_df(get_status_df(ssh_session(),fx_directory()))
    
    status_df_plot<-status_df()
    status_df_plot<-status_df_plot%>%full_join(
      mutate(files_df_fx(), runID=1:nrow(files_df_fx()))%>%
        filter(startsWith(`Engram folder`,"/"))%>%
        select(runID))%>%
      mutate(status=ifelse(is.na(status),"Not processed",status))
    
    ggplot(status_df_plot, aes(status))+
      geom_bar()+theme_bw(18)+ggtitle("Status of sequencing \"runs\"")
    
    
  })
  output$plot2 <- renderPlot({req(QC_df(), status_df()) 
    
    ahmm_files(get_ahmm_files(ssh_session(),fx_directory()))
    
    ahmm_vcf_files_df<-ahmm_files()
    
    ahmm_vcf_files_df<-ahmm_vcf_files_df%>%arrange(mtime)%>%
      mutate(counter = row_number())
    
    ahmm_vcf_files_df<-mutate(ahmm_vcf_files_df, Type="AHMM + VCF generation")
    
    QC_df_plot<-QC_df()
    
    print(head(ahmm_vcf_files_df))
    ahmm_vcf_files_df<-mutate(ahmm_vcf_files_df, mtime=as.POSIXct(as.numeric(mtime), origin = "1970-01-01"))
    print(head(ahmm_vcf_files_df))
    
    QC_df_plot%>%arrange(-Date)%>%distinct(Sample, .keep_all = TRUE)%>%select(Date)%>%
      arrange(Date)%>%mutate(counter = row_number())%>%
      mutate(Type="Initial mapping etc.")%>%
      mutate(Date=as.POSIXct(Date, origin = "1970-01-01"))%>%
      ggplot(aes(Date, counter, color=Type))+
      geom_line()+
      geom_line(data=ahmm_vcf_files_df, aes(mtime, counter, color=Type))+
      theme_bw(18)+ggtitle("Genotyping progress")+
      theme(legend.position="right")+ylab("Number of files")
    
    
    })
  output$plot3 <- renderPlot({ 
    req(ssh_session()) 
    
    today_slurm_df(get_slurm_jobs(ssh_session(),jobs_username()))
    ggplot(today_slurm_df(), aes(jobname, fill=status))+
      geom_bar()+scale_fill_viridis_d()+theme_bw(18)+ggtitle("Today's server work on the Fx pipeline")})
  
  # Handle refresh button
  observeEvent(input$goButton, {
    # Refresh code here
  })
  
  # Show login dialog at start of session
  showModal(modalDialog(
    textInput("username", "Username (you who is running this (needs PW))", value = "user"), # Set default username
    passwordInput("password", "Password", value = ""),
    textInput("jobs_username", "Username (running the pipeline, does not need PW)", value = "user"), # Set default username
    textInput("fx_directory", "FX Directory", value = "/folder/"), # Set default directory
    footer = actionButton("login", "Login"),
    easyClose = FALSE
  ))
}

# Run the application 
shinyApp(ui = ui, server = server)
