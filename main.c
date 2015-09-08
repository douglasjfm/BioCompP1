#include <stdlib.h>
#include <gtk/gtk.h>
#include <gmodule.h>
#include <string.h>

#define UI_FILE "ui.glade"

void alinhar (char *s, char *t)
{
    return;
}

G_MODULE_EXPORT void cb_alinhar (GtkWidget *btn, gpointer data)
{
    char cad1[101], cad2[101], *ret1, *ret2;
    GtkBuilder *builder = GTK_BUILDER(data);
    GtkEntry *uicad1 = NULL, *uicad2 = NULL;
    GtkLabel *label1 = NULL, *label2 = NULL;

    uicad1 = GTK_ENTRY(gtk_builder_get_object(builder,"uiCadeia1"));
    uicad2 = GTK_ENTRY(gtk_builder_get_object(builder,"uiCadeia2"));

    label1 = GTK_LABEL(gtk_builder_get_object(builder,"uiAlignCadeia1"));
    label2 = GTK_LABEL(gtk_builder_get_object(builder,"uiAlignCadeia2"));

    strcpy(cad1,uicad1->text);
    strcpy(cad2,uicad2->text);

    if (!p1_valida_cadeia(cad1) || !p1_valida_cadeia(cad2))
    {
        g_print("Entrada(s) Invalidas.");
    }
    else
    {
        p1_alinhar_s_t(cad1,cad2,&ret1,&ret2);
        gtk_label_set_text(label1, ret1);
        gtk_label_set_text(label2, ret2);
    }
}

int runGTK (int argc, char *argv[])
{
    GtkBuilder *builder = NULL;
    GError *error = NULL;
    GtkWidget *win = NULL;

    /* Initialize GTK+ */
    gtk_init (&argc, &argv);

    builder = gtk_builder_new();

    if (!gtk_builder_add_from_file(builder,UI_FILE,&error))
    {
        g_warning("P1 Erro: %s",error->message);
        g_free(error);
        return 0xE;
    }

    win = GTK_WIDGET(gtk_builder_get_object(builder,"janela"));
    gtk_builder_connect_signals(builder, builder);
    g_signal_connect (win, "destroy", gtk_main_quit, NULL);

    /* Enter the main loop */
    gtk_widget_show_all (win);
    gtk_main ();
    return 0;
}

int main (int argc, char *argv[])
{
    if (argc < 2)
        return runGTK(argc, argv);
    if (argc > 2)
    {
        printf("P1 Console\n");
        if (p1_valida_cadeia(argv[1]) && p1_valida_cadeia(argv[2]))
        {
            printf("Cadeias OK!\n");
            p1_alinhar_s_t(argv[1],argv[2]);
        }
    }
    return 0;
}
